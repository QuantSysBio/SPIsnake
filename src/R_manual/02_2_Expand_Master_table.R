# Define jobs to be executed using Master_table and generated proteome chunks

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master table
Master_table <- fread(file = "Master_table.csv")

# Output dir
directory = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive"

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
  # separate_rows(`Splice type`, sep=",") %>% 
  # mutate(`Splice type`=str_squish(`Splice type`)) %>%
  
  ### N-mers
  mutate(N_mers = str_replace_all(N_mers, pattern = "_", replacement = ",")) %>%
  separate_rows(N_mers, sep=",") %>% 
  mutate(N_mers=as.integer(N_mers)) %>% 
  group_by(Proteome,`Splice type`) %>%
  expand(Proteome, `Splice type`, full_seq(N_mers, 1), Min_Interv_length) %>%
  rename(N_mers = "full_seq(N_mers, 1)",
         Splie_type=`Splice type`) 

Master_table_expanded <- Master_table_expanded %>%
  left_join(proteome_chunks)

fwrite(Master_table_expanded, file = "Master_table_expanded.csv")
