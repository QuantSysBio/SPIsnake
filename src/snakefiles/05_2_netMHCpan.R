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
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
{
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # dir_IC50 = "results/IC50/"
  # cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
  # binders = "results/IC50/IC50_filtered_peptides/PSP_A_11_MeV_BLCL_allFractions_HLA-C07:01.csv.gz"
  # binders = "results/IC50/IC50_filtered_peptides/PSP_L_14_MeV_MA0009-BE08_allFractions_HLA-C07:01.csv.gz"
}

# Experiment_design
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Wildcard
binders <- snakemake@output[[1]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]

### ---------------------------------------------- Predict binding --------------------------------------------------
# CPU
Ncpu = availableCores()
setDTthreads(Ncpu)

# Parameters
Affinity_threshold <- cmd_netMHCpan$Affinity_threshold[str_detect(binders, cmd_netMHCpan$Peptide_file)]

# Run netMHCpan-4.1
cmd <- cmd_netMHCpan$cmds[str_detect(binders, cmd_netMHCpan$Peptide_file)]
print(cmd)
system(cmd, intern = T)
print("Predicted MHC affinity")

### ---------------------------------------------- Select MHC-I binders --------------------------------------------------
netMHCpan_out <- paste0(dir_IC50, "/netMHCpan_output/", cmd_netMHCpan$Peptide_file[str_detect(binders, cmd_netMHCpan$Peptide_file)], ".txt")
print(file.exists(netMHCpan_out))
binders_df <- vroom_lines(file = netMHCpan_out, 
                          num_threads = Ncpu,
                          skip_empty_rows = TRUE,
                          skip = 55)
if (length(binders_df) > 0) {
  binders_df <- binders_df %>%
    str_replace_all(pattern = "[[:space:]]+", " ") %>%
    str_split_fixed(pattern = " ", n = Inf) %>%
    as.data.table() %>%
    lazy_dt() %>%
    rename(MHC = V3, Peptide = V4, `Aff(nM)` = V17) %>%
    select(MHC, Peptide, `Aff(nM)`) %>%
    mutate(`Aff(nM)` = as.numeric(`Aff(nM)`)) %>%
    filter((!is.na(`Aff(nM)`)) & `Aff(nM)` <= Affinity_threshold) %>%
    as.data.table()
  print(head(binders_df))
} else {
  binders_df = data.table(MHC = character(),
                          Peptide = character(),
                          `Aff(nM)` = numeric())
  print("Warning: empty input!")
  print("Make sure that allele is present in netMHCpan-4.1/data/allelenames")
}

IC50_filter_stats <- cmd_netMHCpan[str_detect(binders, cmd_netMHCpan$Peptide_file),] %>%
  mutate(IC50_filtered_peptides = nrow(binders_df))
print(head(IC50_filter_stats))

### ---------------------------------------------- (4) Export --------------------------------------------------
IC50_filter_stats %>%
  vroom_write(append = FALSE, col_names = TRUE, delim = ",",
              file = unlist(snakemake@output[["IC50_filter_stats"]]))

binders_df %>%
  vroom_write(append = FALSE, col_names = TRUE, delim = ",",
              file = unlist(snakemake@output[["binders"]]))
