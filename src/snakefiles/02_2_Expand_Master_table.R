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
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  directory = snakemake@params[["directory"]]
  Master_table <- fread(snakemake@input[["Master_table"]]) %>% as_tibble()
  Ncpu <- snakemake@params[["max_cpus"]]
  
} else {
  ### Manual startup
  Master_table <- fread("Master_table.csv") %>% as_tibble()
  directory = "results/DB_exhaustive"
}

active_Slurm <- ifelse(system("echo $SLURM_JOB_ID") != 0, TRUE, FALSE)
if (active_Slurm) {
  Ncpu <- min(parallelly::availableCores(methods = "Slurm"), Ncpu)
} else {
  Ncpu <- min(parallelly::availableCores())
}

### Find chunks
proteome_chunks <- list.files(paste0(directory, "/Fasta_chunks"), pattern = ".parquet", recursive = FALSE)
proteome_chunks <- lapply(paste0(directory, "/Fasta_chunks/", proteome_chunks), read_parquet, nThread=Ncpu)
proteome_chunks <- lapply(proteome_chunks, setDT)
names(proteome_chunks) <- list.files(paste0(directory, "/Fasta_chunks"), pattern = ".parquet", recursive = FALSE) %>%
  str_remove_all(pattern = ".parquet")
proteome_chunks <- proteome_chunks %>%
  lapply(function(x){
    x %>%
      dplyr::select(tidyr::starts_with("maxE_")) %>%
      as_tibble() %>%
      pivot_longer(names_to = "MaxE", 
                   values_to = "chunk", 
                   cols = tidyr::starts_with("maxE_"), 
                   names_prefix = "maxE_", 
                   names_transform = as.integer) %>%
      unique()
  }) %>%
  rbindlist(idcol = "Proteome") %>%
  as_tibble()

# Tidy format for Master table: PCP & PSP
if (TRUE %in% str_detect(Master_table$enzyme_type, "PCP|PSP|cis-PSP|Cis-PSP")) {
  Master_table_expanded <- Master_table %>% 
    mutate(enzyme_type_lower = str_to_lower(enzyme_type)) %>%
    filter(!str_detect(enzyme_type_lower, str_c(enzymes, collapse = "|"))) %>%
    select(-enzyme_type_lower) %>%
    
    ### N-mers
    mutate(tmp_group = paste(Proteome,`enzyme_type`, N_mers, Max_Interv_length, sep="_")) %>%
    group_by(tmp_group) %>%
    
    mutate(N_mers = str_replace_all(N_mers, pattern = "_", replacement = ",")) %>%
    separate_rows(N_mers, sep=",") %>% 
    mutate(N_mers=as.integer(N_mers)) %>% 
    expand(Proteome, `enzyme_type`, full_seq(N_mers, 1), Max_Interv_length, MaxE, PTMs) %>%
    rename(N_mers = "full_seq(N_mers, 1)") %>%
    mutate(N_mers=as.character(N_mers)) %>% 
    
    ### Splice type
    separate_rows(`enzyme_type`, sep=",") %>% 
    mutate(`enzyme_type`=str_squish(`enzyme_type`)) %>%
    mutate(`enzyme_type`= ifelse(`enzyme_type`=="PSP", 
                                 str_replace_all(`enzyme_type`, pattern = fixed("PSP"), replacement = "cis-PSP"),
                                 `enzyme_type`)) %>%
    mutate(`enzyme_type`=str_replace_all(`enzyme_type`, pattern = "Cis-PSP", replacement = "cis-PSP")) %>%
    
    ### Add peptide chunk files
    left_join(proteome_chunks) %>%
    unique() %>%
    
    ### Create a future wildcard
    mutate(filename = paste(N_mers, enzyme_type, Max_Interv_length, Proteome, chunk, sep = "_")) %>%
    arrange(filename) %>%
    ### Sanity check for redundancies
    ungroup() %>%
    select(-tmp_group) %>%
    unique()
}

# Tidy format for Master table: cleaver
if (TRUE %in% str_detect(str_to_lower(Master_table$enzyme_type), str_c(enzymes, collapse = "|"))) {
  Master_table_enzyme <- Master_table %>%
    mutate(enzyme_type = str_to_lower(enzyme_type)) %>%
    filter(str_detect(enzyme_type, str_c(enzymes, collapse = "|"))) %>%
    
    # Grouping variable
    mutate(tmp_group = paste(Proteome,`enzyme_type`, N_mers, Max_Interv_length, sep="_")) %>%
    group_by(tmp_group) %>%
    
    # Extract enzyme specificity and missed cleavages
    mutate(enzym = str_extract(enzyme_type, pattern = str_c(enzymes, collapse = "|"))) %>%
    mutate(enzyme_type = str_remove_all(enzyme_type, pattern = str_c(enzymes, collapse = "|"))) %>%
    mutate(missedCleavages = as.integer(str_extract(enzyme_type, pattern = "\\d{1,2}"))) %>%
    mutate(missedCleavages = ifelse(is.na(missedCleavages), 0, missedCleavages)) %>%
    mutate(enzyme_type = paste(enzym, missedCleavages, sep = "_")) %>%
    dplyr::select(-c("enzym", "missedCleavages")) %>%
    left_join(proteome_chunks) %>%
    
    ### Create a future wildcard
    mutate(N_mers = as.character(N_mers)) %>% 
    mutate(filename = paste(N_mers, enzyme_type, Max_Interv_length, Proteome, chunk, sep = "_")) %>%
    arrange(filename) %>%
    ### Sanity check for redundancies
    ungroup() %>%
    select(-tmp_group) %>%
    unique()
}

### Common output
if (exists("Master_table_expanded") & exists("Master_table_enzyme")) {
  Master_table_expanded <- bind_rows(Master_table_expanded, 
                                     Master_table_enzyme)
} else if (!exists("Master_table_expanded") & exists("Master_table_enzyme")) {
  Master_table_expanded <- Master_table_enzyme
} 

### Output
if (exists("snakemake")) {
  fwrite(Master_table_expanded, file = unlist(snakemake@output[["Master_table_expanded"]]))
  SPIsnake_log()
  sink()
} else {
  fwrite(Master_table_expanded, file = "results/DB_exhaustive/Master_table_expanded.csv")
}
