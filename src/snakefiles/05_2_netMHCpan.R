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
  # cmd_block_i = "results/IC50/IC50_filtered_peptides/1.csv"
}

# Experiment_design
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Wildcard
cmd_block_i <- snakemake@output[[1]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]

# CPU
Ncpu = availableCores()
setDTthreads(Ncpu)
### ---------------------------------------------- Predict binding --------------------------------------------------
# Select commands
cmd_block_i <- str_split_fixed(cmd_block_i, "/IC50/Seq_stats/", 2)[,2] %>%
  str_remove(".csv")
cmd_netMHCpan <- cmd_netMHCpan %>% 
  filter(cmd_block == cmd_block_i)
cmds <- cmd_netMHCpan %>%
  pull(cmds) %>%
  as.list()

# Run netMHCpan-4.1
mclapply(cmds, mc.cores = Ncpu, function(x){
  print(x)
  system(x, intern = T)
})
print("Predicted MHC affinity")

### ---------------------------------------------- Select MHC-I binders --------------------------------------------------
# Peptide_files <- cmd_netMHCpan$Peptide_file[1:27] %>% as.list()

IC50_filter_stats <- mclapply(X = Peptide_files, mc.cores = Ncpu, FUN = function(x){
  # Parameters
  Affinity_threshold <- cmd_netMHCpan$Affinity_threshold[cmd_netMHCpan$Peptide_file == x]
  
  # Read and filter netMHCpan output
  netMHCpan_out <- paste0(dir_IC50, "/netMHCpan_output/", x, ".txt")
  print(file.exists(netMHCpan_out))
  binders_df <- vroom_lines(file = netMHCpan_out, 
                            num_threads = 1,
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
  # Save filtered peptides
  binders_df %>%
    vroom_write(append = FALSE, col_names = TRUE, delim = ",",
                file = paste0(dir_IC50, "/IC50_filtered_peptides/", x, ".csv.gz"))
  
  # Save stats
  IC50_filter_stats <- cmd_netMHCpan[cmd_netMHCpan$Peptide_file == x,] %>%
    mutate(IC50_filtered_peptides = nrow(binders_df))
  return(IC50_filter_stats)
}) %>%
  bind_rows()
  
### ---------------------------------------------- (4) Export --------------------------------------------------
IC50_filter_stats %>%
  vroom_write(append = FALSE, col_names = TRUE, delim = ",",
              file = unlist(snakemake@output[["IC50_filter_stats"]]))
