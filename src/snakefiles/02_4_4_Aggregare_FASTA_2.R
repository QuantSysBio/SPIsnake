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

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(future))

source("src/snakefiles/functions.R")
cat("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # Master_table_expanded
  Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]]) %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  if (parallelly::supportsMulticore()) {
    plan(multicore)
  } else {
    plan(multisession)
  }
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Folders
  dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- as.logical(snakemake@params[["FASTA_outputs_unfiltered"]])
  FASTA_outputs_MW_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered"]])
  FASTA_outputs_MW_filtered_PTM <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered_PTM"]])
  FASTA_outputs_MW_RT_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_filtered"]])
  FASTA_outputs_MW_RT_IC50_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_IC50_filtered"]])
  
  # Whether to compress CSV output
  compress_CSV <- as.logical(snakemake@params[["compress_CSV"]])
  
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
  if (parallelly::supportsMulticore()) {
    plan(multicore)
  } else {
    plan(multisession)
  }
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Folders
  dir_DB_exhaustive <- "results/DB_exhaustive/"
  dir_DB_out = "results/DB_out/"   
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- TRUE
  FASTA_outputs_MW_filtered <- TRUE
  FASTA_outputs_MW_filtered_PTM <- TRUE
  FASTA_outputs_MW_RT_filtered <- TRUE
  FASTA_outputs_MW_RT_IC50_filtered <- TRUE  
  
  # Whether to compress CSV output
  compress_CSV <- TRUE
  
  # Wildcard
  filename <- "results/DB_out/.Aggregare_FASTA.done"
  filename_out <- filename
}

### Arrow dataset
DB_arrow <- open_dataset(paste0(dir_DB_out, "/arrow/"))
DB_PepProt <- open_dataset(paste0(dir_DB_exhaustive, "/peptide_mapping/")) %>%
  mutate(protein = as.character(protein)) 

### ---------------------------- (2) unfiltered --------------------------------------
if (FASTA_outputs_unfiltered){
  cat("Start saving unfiltered peptide sequences", "\n", as.character(Sys.time()), "\n")
  suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/unfiltered")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/unfiltered/chunks/")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/unfiltered")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/unfiltered/chunks/")))

  groups <- DB_arrow %>%
    select(enzyme, MiSl, proteome) %>%
    unique() %>%
    collect()
  
  for (i in seq_along(groups)) {
    # for (i in 3:4) {
    groups_j <- DB_arrow %>%
      filter(enzyme == groups$enzyme[[i]] & MiSl == groups$MiSl[[i]] & proteome == groups$proteome[[i]]) %>%
      select(length, index) %>%
      unique() %>%
      collect()
    
    if (file.exists(paste0(dir_DB_out, "FASTA/unfiltered/", groups$enzyme[[i]], "_", groups$MiSl[[i]], "_", groups$proteome[[i]], ".fasta"))) {
      file.remove(paste0(dir_DB_out, "FASTA/unfiltered/", groups$enzyme[[i]], "_", groups$MiSl[[i]], "_", groups$proteome[[i]], ".fasta"))
    }
    if (file.exists(paste0(dir_DB_out, "CSV/unfiltered/", groups$enzyme[[i]], "_", groups$MiSl[[i]], "_", groups$proteome[[i]], ifelse(compress_CSV, ".csv.gz", ".csv")))) {
      file.remove(paste0(dir_DB_out, "CSV/unfiltered/", groups$enzyme[[i]], "_", groups$MiSl[[i]], "_", groups$proteome[[i]], ifelse(compress_CSV, ".csv.gz", ".csv")))
    }
    
    for (j in seq_along(groups_j$length)) {
      cat("Starting", groups$enzyme[[i]], " - ", groups$MiSl[[i]],  " - ", groups$proteome[[i]])
      keep_j <- DB_arrow %>%
        filter(enzyme == groups$enzyme[[i]] & MiSl == groups$MiSl[[i]] & proteome == groups$proteome[[i]]) %>%
        filter(length == groups_j$length[[j]] & index == groups_j$index[[j]]) %>%
        select(peptide) %>%
        unique() %>%
        # mutate(fasta = paste0(">pep", "\n", peptide))  %>%
        ungroup() %>%
        # select(fasta) %>%
        collect() %>%
        as.data.table()
      
      if (nrow(keep_j > 0)) {
        # # Save FASTA
        # vroom_write(keep_j, file = paste0(dir_DB_out, "FASTA/unfiltered/", groups$enzyme[[i]], "_", 
        #                                   groups$MiSl[[i]], "_", groups$proteome[[i]], ".fasta"), 
        #             delim = "\n", col_names = F, append = T, quote = "none", num_threads = Ncpu, progress = F)
        
        f <- future({writeXStringSet(AAStringSet(keep_j$peptide), append = T, format = "fasta", 
                          filepath = paste0(dir_DB_out, "FASTA/unfiltered/", groups$enzyme[[i]], "_", 
                                            groups$MiSl[[i]], "_", groups$proteome[[i]], ".fasta"))})
        
        # Save CSV peptide-protein mapping
        DB_PepProt %>%
          filter(enzyme == groups$enzyme[[i]] & MiSl == groups$MiSl[[i]] & proteome == groups$proteome[[i]]) %>%
          filter(length == groups_j$length[[j]] & index == groups_j$index[[j]]) %>%
          select(peptide, protein) %>%
          unique() %>%
          collect() %>%
          as.data.table() %>%
          vroom_write(append = T, delim = ",", num_threads = Ncpu, progress = F,
                      file = paste0(dir_DB_out, "CSV/unfiltered/",
                                    groups$enzyme[[i]], "_", groups$MiSl[[i]], "_", groups$proteome[[i]],
                                    ifelse(compress_CSV, ".csv.gz", ".csv")))
      }
      rm(keep_j)
    }
  }
  cat("Done saving unfiltered peptide sequences", "\n", as.character(Sys.time()), "\n")
}



### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
