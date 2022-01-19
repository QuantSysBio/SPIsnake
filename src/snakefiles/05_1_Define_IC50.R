### ---------------------------------------------- Define IC50  ----------------------------------------------
# description:  Define IC50 prediction and filtering mode
#               
# input:        1. Peptides after MW and RT filtering
#               2. Parameters: alleles and affinities specified in Experiment_design
# output:       
#               - A table with a single line per combination of parameters across calibration datasets
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))
print(sessionInfo())

# {
#   ### setwd("/home/yhorokh/Desktop/wd/tmp/SPIsnake")
#   Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
#   netMHCpan = "bin/netMHCpan-4.1/netMHCpan"
#   dir_DB_PTM_mz = "results/DB_PTM_mz/"
#   dir_IC50 = "results/IC50/"
#   n_netMHCpan_blocks = 3
# }

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (1) Read inputs ----------------------------
netMHCpan <- "/bin/netMHCpan-4.1/netMHCpan"

# Experiment_design
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)
n_netMHCpan_blocks <- snakemake@params[["n_netMHCpan_blocks"]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
{
  suppressWarnings(dir.create(dir_IC50))
  suppressWarnings(dir.create(paste0(dir_IC50, "/netMHCpan_output")))
  suppressWarnings(dir.create(paste0(dir_IC50, "/IC50_filtered_peptides")))
  suppressWarnings(dir.create(paste0(dir_IC50, "/Seq_stats")))
}

### ---------------------------- (2) Define cmds: netMHCpan prediction --------------------------------------
cmds <- tibble(file = list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_for_NetMHCpan"), pattern = ".tsv"),
                   size = file.size(list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_for_NetMHCpan"), pattern = ".tsv", full.names = T))) %>%
  mutate(N_mer = str_split_fixed(file, "_", Inf)[,2]) %>%
  mutate(allele = str_split_fixed(file, pattern = paste0("_", N_mer, "_"), Inf)[,2]) %>%
  mutate(allele = str_split_fixed(allele, pattern = paste0("_ch_"), Inf)[,1]) %>%
  arrange(desc(N_mer), desc(size)) %>% 
  mutate(cmd_block = rep(1:3, length.out = n())) %>%
  arrange(cmd_block, desc(N_mer), desc(size)) %>%
  mutate(cmds = paste(netMHCpan,
                      "-BA", "-inptype 1",
                      "-a", allele,
                      "-l", N_mer, 
                      "-p -f", paste0(dir_DB_PTM_mz, "/unique_peptides_for_NetMHCpan/", file),
                      ">", paste0(dir_IC50,"/netMHCpan_output/", file, ".txt"),
                      "-v"))
t(cmds[1,])

### ---------------------------- (4) Export --------------------------------------
### Peptide_file - column for wildcards
# cmds %>%
#   vroom_write(delim = ",", append = FALSE, col_names = TRUE,
#               file = paste0(dir_IC50, "/cmd_netMHCpan.csv"))

cmds %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file =  unlist(snakemake@output[["cmd_netMHCpan"]]))
