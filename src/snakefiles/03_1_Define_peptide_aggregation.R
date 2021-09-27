### ---------------------------------------------- Define peptide aggregation  ----------------------------------------------
# description:
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               A table with a single line per combination of parameters for peptide across proteome chunks. 
#               They will be used as wildcards to control uniqueness, MW computation, m/z matching and PTM generation. 
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))

{
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake-main")
  # Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv")
  # index_length = 1
  # dir_DB_exhaustive = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_exhaustive"
  # dir_DB_PTM_mz = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_PTM_mz"
}
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read inputs ----------------------------
# Master_table_expanded
Master_table_expanded <- vroom(snakemake@input[["Master_table_expanded"]])

# AA-index
index_length = as.integer(snakemake@params[["AA_index_length"]])
AA = matrix(c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"), 
            nrow=20, ncol = index_length) %>%
  as_tibble() %>% 
  expand.grid() %>%
  as_tibble() %>%
  unite(AA_index, sep = "")

# Output dirs
dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
{
  suppressWarnings(
    dir.create(dir_DB_exhaustive)
  )
  suppressWarnings(
    dir.create(dir_DB_PTM_mz)
  )
  suppressWarnings(
    dir.create(paste0(dir_DB_exhaustive, "/unique_peptide_sequences"))
  )
  suppressWarnings(
    dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched"))
  )
  suppressWarnings(
    dir.create(paste0(dir_DB_PTM_mz, "/files_mz_match"))
  )
  # suppressWarnings(
  #   dir.create(paste0(dir_DB_PTM_mz, "/unique_peptide_MW"))
  # )

  }


### ---------------------------- (2) Define aggregation wildcards --------------------------------------
# Define unique combinations
Peptide_aggregation_table <- Master_table_expanded %>%
  select(N_mers) %>%
  unique() %>%
  expand_grid(AA) %>%
  mutate(AA_length = paste(AA_index, N_mers, sep = "_"))

print(paste("Total", nrow(Peptide_aggregation_table), "combinations of AA_length and N_mers"))

{
  # Check if there are pre-processed proteome chunks
  processed_files <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status"), pattern = ".csv") %>%
    str_remove_all(".csv.gz") %>%
    as_tibble() 
  
  if (nrow(processed_files) > 0) {
   print("Found pre-computed database, adding new chunks")
    
    Process_files <- Master_table_expanded %>%
      filter(!filename %in% processed_files$value) %>%
      left_join(Peptide_aggregation_table) %>%
      unique()
    
    Peptide_aggregation_table <- Master_table_expanded %>%
      filter(!filename %in% processed_files$value) %>%
      left_join(Peptide_aggregation_table) %>%
      select(AA_length) %>%
      unique()
    
    print(paste("Defined aggregation for", nrow(Peptide_aggregation_table), "combinations of AA_length and N_mers"))  
  }
}

### ---------------------------- (3) Export aggregation wildcards --------------------------------------

vroom_write(Peptide_aggregation_table, delim = ",", append = FALSE,
            file = unlist(snakemake@output[["Peptide_aggregation_table"]]))

# vroom_write(Peptide_aggregation_table, delim = ",", append = FALSE,
#             file = paste0(dir_DB_PTM_mz, "/Peptide_aggregation_table.csv"))
