### ---------------------------------------------- Aggregate IC50  ----------------------------------------------
# description:
#               
# input:        1. metadata: Master_table_expanded, cmd_netMHCpan
#               2. Candidate peptide binders filtered from netMHCpan output
#               
# output:       FASTA file per Biological group
#               multimappers 
#                
#               
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))

# {
#   setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
#   Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv", show_col_types = FALSE)
#    Master_table_expanded <- vroom("/home/yhorokh/Downloads/Master_table_expanded.csv")
#   Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
#   Experiment_design$Biological_group <- c(1,2,rep(3,3))
#   dir_DB_exhaustive = "results/DB_exhaustive/"
#   dir_IC50 = "results/IC50/"
#   dir_DB_out = "results/DB_out"
#   cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
#   cmd_netMHCpan <- cmd_netMHCpan %>%
#     mutate(Biological_group = ifelse(Filename == "MeV_BLCL_allFractions", "second", Biological_group)) %>%
#     mutate(Biological_group = ifelse(Filename == "MeV_GRLCL_allFractions", "third", Biological_group))
# }
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (1) Read inputs ----------------------------
# Experiment_design
Master_table_expanded <- vroom(snakemake@input[["Master_table_expanded"]], show_col_types = FALSE)
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Output dirs
dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
dir_DB_out = snakemake@params[["dir_DB_out"]]
dir_IC50 = snakemake@params[["dir_IC50"]]
suppressWarnings(
  dir.create(dir_DB_out)
)

### ---------------------------- (2) Pre-processing --------------------------------------
# Define binder to peptide mapping correspondence
peptide_chunks <- list.files(paste0(dir_IC50, "/IC50_filtered_peptides"), pattern = ".csv.gz", full.names = T) %>%
  as_tibble() %>%
  mutate(Peptide_file = str_remove_all(value, ".csv.gz")) %>%
  mutate(Peptide_file = str_split_fixed(Peptide_file, "/IC50_filtered_peptides/", n = 2)[,2]) %>%
  left_join(cmd_netMHCpan)

peptide_chunks$mapping_prefix <- str_split(peptide_chunks$file, "_", 4, simplify = T) %>%
  as_tibble() %>%
  mutate(V1 = ifelse(str_starts(V1, "PCP"), "PCP_map", "PSP_map")) %>%
  mutate(mapping_prefix = paste(V1, V2, V3, sep="_")) %>%
  pull(mapping_prefix)

# Which predicted binders are non-empty
peptides <- as.list(peptide_chunks$value)
names(peptides) <- peptide_chunks$Peptide_file
if (length(peptides) > 0) {
  peptide_list <- peptides[lapply(peptides, file.size) > 40]
  peptide_chunks <- peptide_chunks %>% filter(Peptide_file %in% names(peptide_list))
  print("Detected non-empty peptide sequences")
}

# Peptide mapping
pep_map <- list.files(paste0(dir_DB_exhaustive, "/peptide_mapping"), pattern = ".csv.gz", full.names = F)
mapping_prefices <- str_c(unique(peptide_chunks$mapping_prefix), collapse = "|")
pep_map <- pep_map[str_starts(pep_map, pattern = mapping_prefices)] %>%
  as.data.table() %>%
  rename(value = ".") %>%
  as.data.table()

pep_map$mapping_prefix <- str_split(pep_map$value, "_", 5, simplify = T) %>%
  as_tibble() %>%
  mutate(mapping_prefix = paste(V1, V2, V3, V4, sep="_")) %>%
  pull(mapping_prefix)

pep_map$pep_length <- str_split_fixed(pep_map$mapping_prefix, "_", 4)[,4] %>% str_length()
pep_map <- split(pep_map, by="mapping_prefix", drop=T, keep.by=T) %>%
  lapply(as_tibble)

### ---------------------------- (3) Informative headers and FASTA export --------------------------------------
mclapply(pep_map, mc.cores = Ncpu, FUN = function(x){
  
  # Identify source proteome and splice type
  df_map <- x %>%
    mutate(Splice_type = ifelse(str_starts(mapping_prefix, "PSP_map_"), "PSP", "PCP")) %>%
    mutate(filename = ifelse(str_ends(value, ".csv.gz"), str_sub(value, end = -nchar(".csv.gz")-1), value)) %>%
    mutate(filename = str_sub(filename, str_length(mapping_prefix) - pep_length + 1))  
  
  mt <- Master_table_expanded %>%
    filter(filename %in% df_map$filename) %>%
    select(filename, Proteome) %>%
    unique()
  df_map <- left_join(df_map, mt)
  
  files_map <- paste0(dir_DB_exhaustive, "/peptide_mapping/", df_map$value) %>% as.list()
  names(files_map) <- paste(df_map$Proteome, df_map$Splice_type, sep = "|")
  
  # Read in peptide map and remove chunk info from the header
  peptide_mapping <- lapply(files_map, vroom, show_col_types = F, num_threads = 1, delim = ",") %>%
    rbindlist(idcol="file")
  peptide_mapping$protein = str_split_fixed(peptide_mapping$protein, "\\|chunk:", 2)[,1]
  peptide_mapping$protein = str_c(peptide_mapping$file, peptide_mapping$protein, sep="|")
  peptide_mapping$file <- NULL
  
  # Predicted binders
  keep_pep <- peptide_chunks$Peptide_file[str_detect(unique(df_map$mapping_prefix), peptide_chunks$mapping_prefix)]
  keep_pep <- peptide_list[names(peptide_list) %in% keep_pep]
  
  # Match predicted binders to mappings for every biological group
  map_biol_group <- peptide_chunks %>%
    filter(Peptide_file %in% names(keep_pep))
  bg <- unique(map_biol_group$Biological_group)
  
  for (biol_group in bg) {
    keep_pep_bg <- keep_pep[names(keep_pep) %in% map_biol_group$Peptide_file[map_biol_group$Biological_group == biol_group]] %>%
      lapply(vroom, show_col_types = F, num_threads = 1, delim = ",") %>%
      rbindlist() %>%
      lazy_dt() %>%
      rename(peptide = Peptide) %>%
      mutate(MW = computeMZ_biostrings(peptide)) %>%
      left_join(peptide_mapping) %>%
      arrange(MW) %>%
      as.data.table() %>%
      unite(col = header, protein, MW, MHC, `Aff(nM)`, remove = T, sep = "|") %>%
      lazy_dt() %>%
      group_by(peptide) %>%
      summarise(header=paste(header,collapse=';')) %>%
      as.data.table()
    
    # Save fasta
    fasta_chunk <- AAStringSet(x = keep_pep_bg$peptide)
    names(fasta_chunk) <- keep_pep_bg$header
    writeXStringSet(fasta_chunk, append = T, format = "fasta",
                    filepath = paste0(dir_DB_out, "/", biol_group, ".fasta"))
  }
})

### ---------------------------- (4) Save stats --------------------------------------
peptide_chunks %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["Summary_stats"]]))
