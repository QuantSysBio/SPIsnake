### ---------------------------------------------- Expand Master table  ----------------------------------------------
# description:  Create a table to control PCP/PSP generation
#               
# input:        1. Master_table contains parameters to be expanded into separate rows
#               2. Assumes that "directory" contains .fasta proteome chunks
#
# output:       
#               - Tidy dataframe with a single combination of proteome chunk and peptide generation parameters per line
#               
# author:       YH

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master table
Master_table <- read.csv(snakemake@input[["Master_table"]])

# Output dir
directory = snakemake@params[["dir_DB_exhaustive"]]

### Find chunks
proteome_chunks <- list.files(directory, pattern = ".fasta", recursive = TRUE) %>%
  strsplit(split = "/", fixed = TRUE) %>%
  as.data.frame() %>%
  t() %>%
  as_tibble() %>%
  rename(Proteome=V1, chunk=V2)

# Tidy format for Master table

Master_table_expanded <- Master_table %>% 
  ### Splice type
  # separate_rows(Splice_type, sep=",") %>% 
  # mutate(Splice_type=str_squish(Splice_type)) %>%
  
  ### N-mers
  mutate(N_mers = str_replace_all(N_mers, pattern = "_", replacement = ",")) %>%
  separate_rows(N_mers, sep=",") %>% 
  mutate(N_mers=as.integer(N_mers)) %>% 
  group_by(Proteome,Splice_type) %>%
  expand(Proteome, Splice_type, full_seq(N_mers, 1), Min_Interv_length) %>%
  rename(N_mers = "full_seq(N_mers, 1)") 

Master_table_expanded <- Master_table_expanded %>%
  left_join(proteome_chunks)


fwrite(Master_table_expanded, file = unlist(snakemake@output[["Master_table_expanded"]]))
