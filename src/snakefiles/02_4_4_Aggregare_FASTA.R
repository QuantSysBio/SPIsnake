### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Create unique .FASTA outputs and their .CSV descriptions
#               
# input:        1. Generate_peptides rules are done
#               2. Peptide sequences are saved as .FASTA and .CSV per chunk. 
#                  For PTM-modified peptides, there should exist arrow datasets
#               3. Experiment_design, Master_table_expanded
# output:       
#               - Unique peptides per Biological group as defined in Experiment_design
#               
# author:       Yehor Horokhovskyi

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # Master_table_expanded
  Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]]) %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Wildcard
  filename <- snakemake@output[[1]]
  filename_out <- filename
  
} else {
  ### Manual startup
  ### setwd("/home/yhorokh/data/Results_reports/wd/SPIsnake-3")
  
  # Master_table_expanded
  Master_table_expanded <- fread("results/DB_exhaustive/Master_table_expanded.csv") %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread("data/Experiment_design.csv") %>% as_tibble()
  
  # CPUs
  Ncpu <- parallelly::availableCores() - 1
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Folders
  dir_DB_out = "results/DB_out/" 
  
  # Wildcard
  filename <- "results/DB_out/.Aggregare_FASTA.done"
  filename_out <- filename
}
### ---------------------------- (2) Aggregate FASTA across chunks --------------------------------------
sub_folders <- c("unfiltered", "MW_filtered", "MW.RT_filtered", "MW.RT.IC50_filtered", "MW_filtered_PTM")
Biological_groups <- unique(Experiment_design$Biological_group)
for (bg in seq_along(Biological_groups)) {
  Biological_group <- Biological_groups[[bg]]
  
  for (i in seq_along(sub_folders)) {
    if (i == 1) {
      FASTA_chunks <- list.files(path = paste0("results/DB_out/FASTA/", sub_folders[[i]], "/chunks/"))
      CSV_chunks <- list.files(path = paste0("results/DB_out/CSV/", sub_folders[[i]], "/chunks/"))
    } else {
      FASTA_chunks <- list.files(path = paste0("results/DB_out/FASTA/", sub_folders[[i]], "/chunks/"), pattern = Biological_group)
      CSV_chunks <- list.files(path = paste0("results/DB_out/CSV/", sub_folders[[i]], "/chunks/"), pattern = Biological_group)
    }
    FASTA_chunks <- FASTA_chunks[str_ends(FASTA_chunks, ".fasta")]
    CSV_chunks <- CSV_chunks[str_ends(CSV_chunks, ".csv|.csv.gz")]
    CSV_chunks_map <- CSV_chunks[str_starts(CSV_chunks, "peptide_mapping_")]
    CSV_chunks <- CSV_chunks[!str_starts(CSV_chunks, "peptide_mapping_")]
    
    ### Concatenate .FASTA
    if (length(FASTA_chunks) > 0) {
      cat("Aggregating FASTA for:", Biological_group, ":", sub_folders[[i]], "\n", as.character(Sys.time()), "\n")
      cmd <- paste0("cat ", str_c(paste0("results/DB_out/FASTA/", sub_folders[[i]], "/chunks/", FASTA_chunks), collapse = " ") ,
                    paste0(" > results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".fasta"))
      system(cmd)
      
      ### Remove duplicates
      cmd <- paste0("seqkit rmdup --by-seq --seq-type protein", 
                    " -j ", Ncpu,
                    " results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".fasta",
                    " -o results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".dup.fasta")
      system(cmd)
      
      ## Check duplicated names
      cmd <- paste0("seqkit rename --by-name",
                    " -j ", Ncpu,
                    " results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".dup.fasta",
                    " -o results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".fasta")
      system(cmd)
      file.remove(paste0("results/DB_out/FASTA/", sub_folders[[i]], "/", Biological_group, ".dup.fasta"))
    }
    
    ### Concatenate .CSV: peptide info
    if (length(CSV_chunks) > 0) {
      cat("Aggregating peptide info .CSV for:", Biological_group, ":", sub_folders[[i]], "\n", as.character(Sys.time()), "\n")
      if (any(str_ends(CSV_chunks, pattern = ".csv.gz"))) {
        CSV_format <- ".csv.gz"
      } else {
        CSV_format <- ".csv"
      }
      cmd <- paste0("csvtk concat ", str_c(paste0("results/DB_out/CSV/", sub_folders[[i]], "/chunks/", CSV_chunks), collapse = " ") ,
                    " -j ", Ncpu, " --keep-unmatched --out-file ",
                    paste0("results/DB_out/CSV/", sub_folders[[i]], "/", Biological_group, CSV_format))
      system(cmd)
    }
    
    ### Concatenate .CSV: peptide-protein mapping
    if (length(CSV_chunks_map) > 0) {
      cat("Aggregating peptide-protein mapping .CSV for:", Biological_group, ":", sub_folders[[i]], "\n", as.character(Sys.time()), "\n")
      if (any(str_ends(CSV_chunks_map, pattern = ".csv.gz"))) {
        CSV_format <- ".csv.gz"
      } else {
        CSV_format <- ".csv"
      }
      cmd <- paste0("csvtk concat ", str_c(paste0("results/DB_out/CSV/", sub_folders[[i]], "/chunks/", CSV_chunks_map), collapse = " ") ,
                    " -j ", Ncpu, " --keep-unmatched --out-file ",
                    paste0("results/DB_out/CSV/", sub_folders[[i]], "/peptide_mapping_", Biological_group, CSV_format))
      system(cmd)
    }
  } # End sub-folders
} # End Biological group

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
