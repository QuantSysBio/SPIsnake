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


### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
print(sessionInfo())

### ---------------------------- (1) Read input file and extract info ----------------------------
# {
#   ### Manual startup
#   Master_table <- read.csv("Master_table.csv")
#   directory = "/home/yhorokh/SNAKEMAKE/QSB_SPIsnake/results/DB_exhaustive"
# }

# Master table
Master_table <- read.csv(snakemake@input[["Master_table"]])

# Output dir
directory = snakemake@params[["directory"]]
print(directory)

### Find chunks
proteome_chunks <- list.files(paste0(directory, "/Fasta_chunks"), pattern = ".fasta", recursive = TRUE) %>%
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
    mutate(tmp_group = paste(Proteome,`Splice_type`, N_mers, sep="_")) %>%
    group_by(tmp_group) %>%
    
    mutate(N_mers = str_replace_all(N_mers, pattern = "_", replacement = ",")) %>%
    separate_rows(N_mers, sep=",") %>% 
    mutate(N_mers=as.integer(N_mers)) %>% 
    expand(Proteome, `Splice_type`, full_seq(N_mers, 1), Min_Interv_length, MaxE) %>%
    rename(N_mers = "full_seq(N_mers, 1)") %>%
    
    ### Splice type
    separate_rows(`Splice_type`, sep=",") %>% 
    mutate(`Splice_type`=str_squish(`Splice_type`)) %>%
    mutate(`Splice_type`= ifelse(`Splice_type`=="PSP", 
            str_replace_all(`Splice_type`, pattern = fixed("PSP"), replacement = "cis-PSP"),
            `Splice_type`)) %>%
    mutate(`Splice_type`=str_replace_all(`Splice_type`, pattern = "Cis-PSP", replacement = "cis-PSP")) %>%
    
    ### Attributes to keep from Master_table and peptide chunk files
    left_join(select(Master_table, Proteome, PTMs)) %>%
    left_join(proteome_chunks) %>%
    
    ### Create a future wildcard
    mutate(filename = paste(N_mers, Splice_type, Min_Interv_length, chunk, sep = "_")) %>%
    arrange(filename) %>%
    ### Sanity check for redundancies
    ungroup() %>%
    select(-tmp_group) %>%
    unique() 
  }

# Tidy format for indices
### PSP
PSP_indices <- Master_table_expanded %>% 
  filter(`Splice_type` == "cis-PSP") %>% 
  ungroup() %>%
  
  # Create a future wildcard
  mutate(PSP_index = paste(N_mers, Min_Interv_length, sep = "_")) %>%
  select(PSP_index) %>%
  unique() 

### Output
fwrite(Master_table_expanded, file = unlist(snakemake@output[["Master_table_expanded"]]))
fwrite(PSP_indices, file = unlist(snakemake@output[["PSP_indices"]]))
