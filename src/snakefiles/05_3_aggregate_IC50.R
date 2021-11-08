### ---------------------------------------------- Define RT train  ----------------------------------------------
# description:
#               
# input:        1. RT calibration datasets: sequence + time (seconds)
#               2. Parameters: n_folds, train_proportion
# output:       
#               A table with a single line per combination of parameters across calibration datasets. 
#               Data split into train/test sets
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))

{
  setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
  dir_IC50 = "results/IC50/"
  dir_DB_out = "results/DB_out"
  cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
  informative_headers = TRUE
}
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (1) Read inputs ----------------------------
# Experiment_design
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Output dirs
dir_DB_out = snakemake@params[["dir_DB_out"]]
dir_IC50 = snakemake@params[["dir_IC50"]]
suppressWarnings(
  dir.create(dir_DB_out)
)

### ---------------------------- (2) Save fasta --------------------------------------
# Select all files with the same {AA-index}{length}
peptide_chunks <- list.files(paste0(dir_IC50, "/IC50_filtered_peptides"), pattern = ".csv.gz", full.names = T) %>%
  as_tibble() %>%
  mutate(Peptide_file = str_remove_all(value, ".csv.gz")) %>%
  mutate(Peptide_file = str_split_fixed(Peptide_file, "/IC50_filtered_peptides/", n = 2)[,2]) %>%
  left_join(cmd_netMHCpan)

peptides <- as.list(peptide_chunks$value)
names(peptides) <- peptide_chunks$Peptide_file

if (length(peptides) > 0) {
  peptide_list <- peptides %>%
    mclapply(FUN = vroom, delim = ",", mc.cores = Ncpu, show_col_types = FALSE) 
  print("Reading in peptide sequences: Done")
}

### ---------------------------- (2) Informative headers --------------------------------------
if (informative_headers == TRUE) {
  
}


DB


### ---------------------------- (4) Save outputs --------------------------------------
Summary_stats %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["Summary_stats"]]))