### ---------------------------------------------- Predict_MHC-I_affinity ----------------------------------------------
# description:  Predict MHC-I affinity for a given combination of alleles and save the filtered peptides 
#               
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .csv.gz
#               - Unique peptides MW .csv.gz
#               - Unique peptides after m/z and RT matching per experiment .csv.gz
#               - Peptide-mass_list matching
#               - Peptide filtering stats .csv
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
{
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # dir_IC50 = "results/IC50/"
  # cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
  # 
  # binders = "results/IC50/IC50_filtered_peptides/PCP_T_15_MeV_MA0009-BE08_allFractions_H-2-Kb.csv.gz"
}

# Experiment_design
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Wildcard
binders <- snakemake@output[[1]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]

### ---------------------------------------------- Predict binding --------------------------------------------------
# Parameters
Affinity_threshold <- cmd_netMHCpan$Affinity_threshold[str_detect(binders, cmd_netMHCpan$Peptide_file)]

# Run netMHCpan-4.1
cmd <- cmd_netMHCpan$cmds[str_detect(binders, cmd_netMHCpan$Peptide_file)]
print(cmd)
system(cmd, intern = T)
print("Predicted MHC affinity")

suppressWarnings(
  netMHCpan_out <- read_delim(
    file = paste0(dir_IC50, "/netMHCpan_output/", cmd_netMHCpan$Peptide_file[str_detect(binders, cmd_netMHCpan$Peptide_file)], ".txt"),
    trim_ws = TRUE,
    delim = " ", 
    quote = "#", 
    skip = 55, 
    col_names = c("Pos", "MHC", "Peptide", "Core", 
                  "Of", "Gp", "Gl", "Ip", "Il", "Icore", "Identity", 
                  "Score_EL", "Perc_Rank", "Score_BA", "Perc_Rank_BA", "Aff(nM)", "BindLevel2", "BindLevel"),
    show_col_types = FALSE) %>%
    select(MHC, Peptide, `Aff(nM)`) %>%
    mutate(`Aff(nM)` = as.numeric(`Aff(nM)`)) %>%
    filter(!MHC %in% c("PEPLIST","PEPLIST.")) %>%
    filter(!is.na(MHC))
)
print(head(netMHCpan_out))

### ---------------------------------------------- Select MHC-I binders --------------------------------------------------
binders_df <- netMHCpan_out %>%
  select(-MHC) %>%
  filter(`Aff(nM)` <= Affinity_threshold)
print(head(binders_df))

IC50_filter_stats <- cmd_netMHCpan[str_detect(binders, cmd_netMHCpan$Peptide_file),] %>%
  mutate(IC50_filtered_peptides = nrow(binders_df))
print(head(IC50_filter_stats))

### ---------------------------------------------- (4) Export --------------------------------------------------
IC50_filter_stats %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["IC50_filter_stats"]]))

binders_df %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["binders"]]))

