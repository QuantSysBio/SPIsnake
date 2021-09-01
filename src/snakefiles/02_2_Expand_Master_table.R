# Define jobs to be executed using Master_table and generated proteome chunks

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master table
Master_table <- read.csv(snakemake@input[["Master_table"]])
print(Master_table)

# Output dir
directory = snakemake@params[["directory"]]
print(directory)


### Find chunks
proteome_chunks <- list.files(directory, pattern = ".fasta", recursive = TRUE) %>%
  strsplit(split = "/", fixed = TRUE) %>%
  as.data.frame() %>%
  t() %>%
  as_tibble() %>%
  rename(Proteome=V1, chunk=V2) %>%
  mutate(MaxE = str_split_fixed(chunk, pattern = fixed("_"), n = 2)[,1]) %>%
  mutate(chunk = str_split_fixed(chunk, pattern = fixed("_"), n = 2)[,2]) %>%
  mutate_at("MaxE", as.numeric) 

# Tidy format for Master table
{
  Master_table_expanded <- Master_table %>% 
    ### N-mers
    mutate(N_mers = str_replace_all(N_mers, pattern = "_", replacement = ",")) %>%
    separate_rows(N_mers, sep=",") %>% 
    mutate(N_mers=as.integer(N_mers)) %>% 
    group_by(Proteome,`Splice_type`) %>%
    expand(Proteome, `Splice_type`, full_seq(N_mers, 1), Min_Interv_length, MaxE) %>%
    rename(N_mers = "full_seq(N_mers, 1)") %>%
    
    ### Splice type
    separate_rows(`Splice_type`, sep=",") %>% 
    mutate(`Splice_type`=str_squish(`Splice_type`)) %>%
    
    ### Attributes to keep from Master_table and peptide chunk files
    left_join(select(Master_table, Proteome, PTMs)) %>%
    left_join(proteome_chunks) %>%
    ### Create a future wildcard
    mutate(filename = paste(Splice_type, N_mers, PTMs, Min_Interv_length, chunk, sep = "_")) %>%
    arrange(filename) %>%
    ### Sanity check for redundancies
    unique() 
  }


# Output
print(head(Master_table_expanded))
fwrite(Master_table_expanded, file = unlist(snakemake@output[["Master_table_expanded"]]))
