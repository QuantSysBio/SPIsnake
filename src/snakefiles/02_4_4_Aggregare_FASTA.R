### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Create unique .FASTA outputs
#               
# input:        1. Generate_peptides rules are done
#               2. Peptide sequences are saved as .FASTA and .CSV per chunk. 
#                  For PTM-modified peptides, there should exist arrow datasets
#               3. Experiment_design, Master_table_expanded
# output:       
#               - Unique peptides per Biological group as defined in Experiment_design
#               
# author:       YH

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # Master_table_expanded
  Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]]) %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- as.logical(snakemake@params[["FASTA_outputs_unfiltered"]])
  FASTA_outputs_MW_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered"]])
  FASTA_outputs_MW_filtered_PTM <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered_PTM"]])
  FASTA_outputs_MW_RT_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_filtered"]])
  FASTA_outputs_MW_RT_IC50_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_IC50_filtered"]])
  
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Folders
  dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  duckdb_RAM = snakemake@params[["duckdb_RAM"]]
  duckdb_max_retries = snakemake@params[["duckdb_max_retries"]]
  duckdb_max_filesize = snakemake@params[["duckdb_max_filesize"]]
  
  # Wildcard
  filename <- snakemake@output[[1]]
} else {
  ### Manual startup
  
  # Master_table_expanded
  Master_table_expanded <- fread("results/DB_exhaustive/Master_table_expanded.csv") %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread("data/Experiment_design.csv") %>% as_tibble()
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- TRUE
  FASTA_outputs_MW_filtered <- TRUE
  FASTA_outputs_MW_filtered_PTM <- TRUE
  FASTA_outputs_MW_RT_filtered <- TRUE
  FASTA_outputs_MW_RT_IC50_filtered <- TRUE  
  
  # CPUs
  Ncpu <- parallelly::availableCores() - 1
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  ### duckdb settings
  duckdb_RAM <- 0.8
  duckdb_max_retries = 20
  duckdb_max_filesize = 2
  
  # Folders
  dir_DB_out = "results/DB_out/" 
  
  # Wildcard
  filename <- "results/DB_out/.Aggregare_FASTA.done"
  filename_out <- filename
}

##№ duckdb parameters
duckdb_max_filesize <- duckdb_max_filesize * (2^30)
cat(as.character(Sys.time()), " - ", "duckdb max filesize: ", duckdb_max_filesize, "\n")

aggregation_batch_i <- str_remove(filename, ".done")
aggregation_batch_i <- str_split_fixed(aggregation_batch_i, "results/DB_out/.", 2)[,2][[1]]

duckdb_temp_dir <- paste0(dir_DB_out, "/duckdb/tmp/unique_peptides_", aggregation_batch_i, "/")
table_name = paste0(dir_DB_out, "/duckdb/databases/", "duckdb_aggregate_unique_", aggregation_batch_i)

Max_RAM <- system("free -g", intern = T)[[2]] %>% 
  gsub(pattern = "\\s+", replacement = " ") %>% 
  str_split_i(pattern = " ", i = 2) %>% 
  as.numeric()

timeout <- 90
duckdb_RAM <- Max_RAM * duckdb_RAM

##№ ----------------------------------- Define groups -----------------------------------
### Arrow dataset
DB_arrow <- open_dataset(paste0(dir_DB_out, "/arrow/"))

DB_arrow_files <- DB_arrow$files %>%
  as_tibble() %>%
  mutate(filename = value) %>%
  mutate(file_size = file.size(filename)) %>%
  mutate(value = str_split_fixed(value, "/arrow/", 2)[,2]) %>%
  mutate(value = str_split_fixed(value, "/part-", 2)[,1]) %>%
  mutate(value = str_remove_all(value, "enzyme=|MiSl=|proteome=|index=|length=")) %>%
  separate(value, into = c("index", "length", "proteome", "enzyme", "MiSl"), sep = "/")  %>%
  mutate(MiSl = as.integer(MiSl),
         length = as.integer(length))
DB_arrow_files

groups_j <- DB_arrow_files %>%
  select(index, length) %>%
  unique() %>%
  arrange(index, length)
groups_j

groups <- DB_arrow_files %>%
  select(enzyme, MiSl, proteome) %>%
  unique()
groups

### ---------------------------- (2) unfiltered --------------------------------------
FASTA_outputs_unfiltered <- FALSE
if (FASTA_outputs_unfiltered){
  cat(as.character(Sys.time()), " - ", "Start saving unfiltered peptide sequences", "\n")
  dir.create(paste0(dir_DB_out, "/FASTA/unfiltered"), recursive = T, showWarnings = T)
  file.remove(paste0(dir_DB_out, "FASTA/.unfiltered.done"))
  
  make_fasta_batched(DB_arrow_files = DB_arrow_files,
                     groups = groups,
                     groups_j = groups_j,
                     timeout = 10,
                     retries = 100,
                     wait_on_restart = 0,
                     filter_prefix = NULL, 
                     filter_suffix = NULL,
                     fasta_Biological_group = "",
                     path_out = paste0(dir_DB_out, "FASTA/unfiltered/unfiltered"))
  
  file.create(paste0(dir_DB_out, "FASTA/.unfiltered.done"))
  cat(as.character(Sys.time()), " - ", "Done saving unfiltered peptide sequences","\n")
} # End FASTA unfiltered

### ---------------------------- (3) Filtered FASTA --------------------------------------
# Experiment_design$Biological_group <- c("K562")
# Experiment_design$Biological_group <- c("R1", "R2")

Biological_groups <- unique(Experiment_design$Biological_group)
for (bg in seq_along(Biological_groups)) {
  Biological_group_bg <- Biological_groups[[bg]] 
  select_mass_list <- Experiment_design %>%
    filter(Biological_group == Biological_group_bg) %>%
    as_tibble() %>%
    pull(Filename)
  print(Biological_group_bg)
  print(select_mass_list)
  
  if (FASTA_outputs_MW_filtered){
    if (any(str_starts(names(open_dataset(DB_arrow_files$filename[[1]])), "MW.exists"))) {
      cat(as.character(Sys.time()), " - ", "Start saving MW filtered peptide sequences", "\n")
      dir.create(paste0(dir_DB_out, "/FASTA/MW_filtered"), recursive = T, showWarnings = T)
      file.remove(paste0(dir_DB_out, "FASTA/.MW_filtered.done"))
      
      make_fasta_batched(DB_arrow_files = DB_arrow_files,
                         groups = groups,
                         groups_j = groups_j,
                         timeout = 7,
                         retries = 100,
                         wait_on_restart = 0,
                         filter_prefix = "MW.exists:", 
                         filter_suffix = select_mass_list,
                         fasta_Biological_group = Biological_groups[[bg]],
                         path_out = paste0(dir_DB_out, "FASTA/MW_filtered/MW_filtered"))
      
      file.create(paste0(dir_DB_out, "FASTA/.MW_filtered.done"))
      cat(as.character(Sys.time()), " - ", "Done saving MW filtered peptide sequences","\n")
    }
  }
  
  if (FASTA_outputs_MW_RT_filtered){
    if (any(str_starts(names(open_dataset(DB_arrow_files$filename[[1]])), "MW.RT.exists"))) {
      cat(as.character(Sys.time()), " - ", "Start saving MW-RT filtered peptide sequences", "\n")
      dir.create(paste0(dir_DB_out, "/FASTA/MW_RT_filtered"), showWarnings = T, recursive = T)
      file.remove(paste0(dir_DB_out, "FASTA/.MW_RT_filtered.done"))
      
      make_fasta_batched(DB_arrow_files = DB_arrow_files,
                         groups = groups,
                         groups_j = groups_j,
                         timeout = 10,
                         retries = 100,
                         wait_on_restart = 0,
                         filter_prefix = "MW.RT.exists:", 
                         filter_suffix = select_mass_list,
                         fasta_Biological_group = Biological_groups[[bg]],
                         path_out = paste0(dir_DB_out, "FASTA/MW_RT_filtered/MW_RT_filtered"))
      
      file.create(paste0(dir_DB_out, "FASTA/.MW_RT_filtered.done"))
      cat(as.character(Sys.time()), " - ", "Done saving MW-RT filtered filtered peptide sequences","\n")
    }
  }
  
  if (FASTA_outputs_MW_RT_IC50_filtered){
    if (any(str_starts(names(open_dataset(DB_arrow_files$filename[[1]])), "Predicted_binder:"))) {
      cat(as.character(Sys.time()), " - ", "Start saving MW-RT-IC50 filtered peptide sequences", "\n")
      dir.create(paste0(dir_DB_out, "/FASTA/MW_RT_IC50"), showWarnings = T, recursive = T)
      file.remove(paste0(dir_DB_out, "FASTA/.MW_RT_IC50_filtered.done"))
      
      make_fasta_batched(DB_arrow_files = DB_arrow_files,
                         groups = groups,
                         groups_j = groups_j,
                         timeout = 10,
                         retries = 100,
                         wait_on_restart = 0,
                         filter_prefix = "Predicted_binder:", 
                         filter_suffix = select_mass_list,
                         fasta_Biological_group = Biological_groups[[bg]],
                         path_out = paste0(dir_DB_out, "FASTA/MW_RT_IC50/MW_RT_IC50"))
      
      file.create(paste0(dir_DB_out, "FASTA/.MW_RT_IC50.done"))
      cat(as.character(Sys.time()), " - ", "Done saving MW-RT-IC50 filtered filtered peptide sequences","\n")
    }
  }
  if (FASTA_outputs_MW_filtered_PTM){
    # TODO implement
  }
} # End Biological group

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
