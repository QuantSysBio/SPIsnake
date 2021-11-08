### ---------------------------------------------- Define IC50  ----------------------------------------------
# description:  Define IC50 prediction and filtering mode
#               
# input:        1. Peptides after MW and RT filtering
#               2. Parameters: alleles and affinities specified in Experiment_design
# output:       
#               A table with a single line per combination of parameters across calibration datasets
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))

# {
#   setwd("/home/yhorokh/Snakemake/SPIsnake/")
#   Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
#   netMHCpan = "/bin/netMHCpan-4.1/netMHCpan"
#   dir_DB_PTM_mz = "results/DB_PTM_mz/"
#   dir_IC50 = "results/IC50/"
# }

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
netMHCpan <- snakemake@input[["netMHCpan"]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]

# RT-filtered peptides
peptides <- tibble(file = list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched"))) %>%
  mutate(N_mer = str_split_fixed(file, "_", Inf)[,3]) %>%
  mutate(Filename = str_remove_all(str_split_fixed(file, "_", 4)[,4], pattern = ".tsv"))
peptides

{
  suppressWarnings(
    dir.create(dir_IC50)
  )
  suppressWarnings(
    dir.create(paste0(dir_IC50, "/netMHCpan_output"))
  )
  suppressWarnings(
    dir.create(paste0(dir_IC50, "/IC50_filtered_peptides"))
  )
  suppressWarnings(
    dir.create(paste0(dir_IC50, "/Seq_stats"))
  )
}

### ---------------------------- (2) Define IC50 prediction --------------------------------------
IC50_aggregation_table <- Experiment_design %>%
  # tidy format
  tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]") %>%
  mutate(`MHC-I_alleles` = str_squish(`MHC-I_alleles`),
         Affinity_threshold = str_squish(Affinity_threshold)) 
IC50_aggregation_table

### ---------------------------- (3) Define cmds: netMHCpan prediction --------------------------------------
cmds <- IC50_aggregation_table %>%
  left_join(peptides) %>%
  mutate(Peptide_file = str_replace_all(file, ".tsv", paste0("_",`MHC-I_alleles`)))
head(cmds$Peptide_file)

cmds$cmds <- paste(netMHCpan,
                   "-BA", "-inptype 1",
                   "-a", cmds$`MHC-I_alleles`,
                   "-l", cmds$N_mer, 
                   "-p -f", paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", cmds$file),
                   ">", paste0("results/IC50/netMHCpan_output/", cmds$Peptide_file, ".txt"),
                   "-v")

### ---------------------------- (4) Export --------------------------------------
# Peptide_file - column for wildcards
{
  IC50_aggregation_table  %>%
    vroom_write(delim = ",", append = FALSE, col_names = TRUE,
                file = paste0(dir_IC50, "/IC50_aggregation_table.csv"))
  
  cmds %>%
    vroom_write(delim = ",", append = FALSE, col_names = TRUE,
                file = paste0(dir_IC50, "/cmd_netMHCpan.csv"))
}

IC50_aggregation_table  %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["IC50_aggregation_table"]]))

cmds %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file =  unlist(snakemake@output[["cmd_netMHCpan"]]))

