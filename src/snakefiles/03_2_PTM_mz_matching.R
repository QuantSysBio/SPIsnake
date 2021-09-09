### ---------------------------------------------- Aggregate_PTM_mz_matching ----------------------------------------------
# description:  Find a unique set of peptide sequences, compute molecular weight (MW) and do m/z matching with all input mass_lists. 
#               If PTMs are required, generate them too. 
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .csv.gz
#               - Unique peptides MW .csv.gz
#               - Unique peptides after m/z matching per experiment .csv.gz
#               - Peptide-mass_list matching
#               - Peptide filtering stats .csv
#               
# author:       YH, JL, KP

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(Biostrings)
library(data.table)
library(dtplyr)
library(dplyr)
library(seqinr)
library(stringr)
library(parallel)
library(parallelly)
library(vroom)
source("src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master_table_expanded
Master_table_expanded <- vroom(snakemake@input[["Master_table_expanded"]])
# Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv")

# Master_table_expanded
Peptide_aggregation_table <- vroom(snakemake@input[["Peptide_aggregation_table"]], delim = ",")
# Peptide_aggregation_table <- vroom("results/DB_PTM_mz/Peptide_aggregation_table.csv", delim = ",")

# Output dir
{
  dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
# dir_DB_exhaustive = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_exhaustive"

dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
# dir_DB_PTM_mz = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_PTM_mz"
}


### Wildcard
# filename = "results/DB_PTM_mz/chunk_aggregation_status/S_8.csv"
filename = snakemake@output[[1]]
filename <- filename %>%
  str_split_fixed(pattern = fixed("chunk_aggregation_status/"), n = 2)
filename <- filename[,2] %>%
  str_remove(pattern = ".csv")
print(filename)


### Mass_lists
MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
  as_tibble() %>%
  mutate(file = str_remove_all(value, ".txt")) 

PCP <- peptides[str_detect(peptides, "/PCP_")]  %>%
  lapply(FUN = vroom, delim = ",", num_threads = Ncpu, show_col_types = FALSE) %>%
  rbindlist(idcol = "file") %>%
  lazy_dt()

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

# Save into chunks according to first N letters
index_length = as.integer(snakemake@params[["AA_index_length"]])
# index_length = 1

### ---------------------------- (2) Operation mode --------------------------------------
# Check if there exist previous outputs to be updated:
processed_files <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status"), pattern = ".csv") %>%
  str_remove_all(".csv.gz") %>%
  as_tibble() 

operation_mode = ifelse(nrow(processed_files) > 0, "Update", "Generation")

### ---------------------------- (3) Uniqueness --------------------------------------
{
  # Select all files with the same {AA-index}{length}
  peptide_chunks <- list.files(paste0(dir_DB_exhaustive, "/peptide_seqences"), pattern = ".csv.gz", full.names = T) %>%
    as_tibble() %>%
    mutate(file = str_remove_all(value, ".csv.gz")) %>%
    mutate(file = str_split_fixed(file, "/peptide_seqences/", n = 2)[,2]) %>%
    mutate(Splice_type = ifelse(str_starts(file, "PCP_"), "PCP", "other")) %>%
    mutate(Splice_type = ifelse(str_starts(file, "PSP_"), "PSP", Splice_type)) %>%
    mutate(file = str_sub(file, start=5)) %>%
    mutate(file = str_remove(file, "PCP_")) %>%
    mutate(file = str_remove(file, "cis-PSP_")) %>%
    filter(str_starts(file, filename))
  
  peptides <- as.list(peptide_chunks$value)
  names(peptides) <- peptide_chunks$file
  
  # PCP
  PCP <- peptides[str_detect(peptides, "/PCP_")]  %>%
    lapply(FUN = vroom, delim = ",", num_threads = Ncpu, show_col_types = FALSE) %>%
    rbindlist(idcol = "file") %>%
    lazy_dt()
  
  # PSP
  PSP <- peptides[str_detect(peptides, "/PSP_")]  %>%
    lapply(FUN = vroom, delim = ",", num_threads = Ncpu, show_col_types = FALSE) %>%
    rbindlist(idcol = "file") %>%
    lazy_dt()
}

### ---------------------------- (4) Compute MW --------------------------------------
### Make a unique set of peptides
{
  PCP <- PCP %>%
    select(peptide) %>%
    unique() %>%
    mutate(MW=computeMZ_biostrings(peptide)) %>%
    as.data.table()
}
{
  PSP <- PSP %>%
    select(peptide) %>%
    unique() %>%
    mutate(MW=computeMZ_biostrings(peptide)) %>%
    as.data.table()
}

### ---------------------------- (5) Generate PTMs --------------------------------------




### ---------------------------- (6) m/z matching --------------------------------------
# Select !none
input <- PCP

for (i in 1:nrow(MS_mass_lists)) {
  print(MS_mass_lists$value[i])
  
  tolerance = 20
  mzList = vroom(file = paste0("data/MS_mass_lists/", MS_mass_lists$value[i]), delim = ",", col_names = "Precursor_mass") %>%
    lazy_dt() %>%
    mutate(Min = Precursor_mass - Precursor_mass * tolerance * 10 ** (-6)) %>%
    mutate(Max = Precursor_mass + Precursor_mass * tolerance * 10 ** (-6)) %>%
    select(-Precursor_mass) %>%
    as.data.table()
  
  # Which peptide sequences pass the MW filter
  index_subset <- mixANDmatch3(mzMin=mzList$Min, mzMax=mzList$Max, MW0 = PCP$MW)
  input <- input[index_subset,]
}




### ---------------------------- (7) Export --------------------------------------





