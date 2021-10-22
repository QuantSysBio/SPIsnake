### ---------------------------------------------- Generate PCP/PSP ----------------------------------------------
# description:  Generate unique PCP/PSP sequences
#               
# input:        1. Confirmed existance of Proteome_chunk.fasta
#               2. Parameters: index_length, peptide length, min_intervening_sequence length, output directory
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences .csv.gz
#               - Protein-peptide mapping .csv.gz
#               - Peptide generation stats .csv
#               
# author:       YH, JL, KP

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

# Functions
seq_list_to_dt <- function(seq_list){
  rbindlist(lapply(seq_list, as.data.table), idcol = "id")
}

# Manual startup
{
  # ##setwd("/home/yhorokh/Snakemake/SPIsnake-main")
  # directory = "results/DB_exhaustive/"
  # dir_DB_Fasta_chunks = "results/DB_exhaustive/Fasta_chunks/"
  # Master_table_expanded <- read.csv("results/DB_exhaustive/Master_table_expanded.csv")
  # filename = "results/DB_exhaustive/Seq_stats/8_cis-PSP_Unimod_25_Measles_CDS_6_frame_1_48.fasta.csv.gz"
  # index_length = 1
  # max_protein_length = 100
  # 
  # 
  # filename <- str_remove(filename, ".csv.gz") %>%
  #   str_split_fixed(pattern = fixed("Seq_stats/"), n = 2)
  # filename <- filename[,2]
  # print(filename)
  # 
  # params <- Master_table_expanded[Master_table_expanded$filename == filename,]
  # print(t(params))
}

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master_table_expanded
Master_table_expanded <- read.csv(snakemake@input[["Master_table_expanded"]])

# Fasta input
dir_DB_Fasta_chunks = snakemake@params[["dir_DB_Fasta_chunks"]]

# Output dir
directory = snakemake@params[["directory"]]
suppressWarnings(
  dir.create(paste0(directory, "/peptide_seqences"))
)
suppressWarnings(
  dir.create(paste0(directory, "/peptide_mapping"))
)

# Wildcard
filename = snakemake@output[[1]]
{
  filename <- str_remove(filename, ".csv.gz") %>%
    str_split_fixed(pattern = fixed("Seq_stats/"), n = 2)
  filename <- filename[,2]
  print(filename)
}

params <- Master_table_expanded[Master_table_expanded$filename == filename,]
print(t(params))

### Extract parameters
# Input fasta
proteome = list.files(dir_DB_Fasta_chunks, pattern = params$chunk, recursive = T, full.names = T)
proteome = proteome[str_ends(proteome, ".fasta")]
dat = read.fasta(file=proteome, 
                 seqtype="AA", as.string = TRUE)

# Keep only proteome name
proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]

# Nmers
Nmers = as.numeric(params$N_mers)

# max intervening sequence length
MiSl = as.numeric(params$Min_Interv_length)

# Splice type
Splice_type = as.character(params$Splice_type)

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
data.table::setDTthreads(Ncpu)

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X|\\*)"

# Save into chunks according to first N letters
index_length = as.integer(snakemake@params[["AA_index_length"]])

# Save into chunks according to first N letters
max_protein_length = as.integer(snakemake@params[["max_protein_length"]])

# Load pre-computed index to speed-up PSP generation
index_list_result = paste(Nmers, MiSl, sep = "_")
if (grepl("cis-PSP", Splice_type)==T) {
  load(paste0(directory, "/PSP_indices/", index_list_result, ".rds"))
}


### ---------------------------- (2) Compute PCP and PSP --------------------------------------
print(Sys.time())
print(paste("Starting length: ", Nmers))

# Compute PCP and PSP per sequence
dat_sort <- dat[names(sort(unlist(lapply(dat, nchar)), decreasing = T))]

if (grepl("cis-PSP", Splice_type)==T) {
  print("Computing PCP and cis-PSP")
  
  Pep_list <- mclapply(dat_sort, 
                       mc.cores = Ncpu,
                       CutAndPaste_seq_from_big_sp_fast, 
                       big_sp_input=index_list_result,
                       nmer = Nmers, 
                       MiSl=MiSl)
  print(Sys.time())
  print("Computed PCP/PSP")
} else if (Splice_type == "PCP") {
  print("Computing PCP only")
  
  Pep_list <- mclapply(dat_sort, mc.cores = Ncpu, CutAndPaste_seq_PCP, nmer = Nmers)
  print(Sys.time())
  print("Computed PCP")
} else if (!exists("Pep_list")) {
  print("Unknown splicing type, no peptides will be generated")
} 

# Tidy format
PSP <- seq_list_to_dt(lapply(Pep_list, `[[`, 1))
PCP <- seq_list_to_dt(lapply(Pep_list, `[[`, 2))
Seq_stats <- rbindlist(lapply(Pep_list, `[[`, 3))
print("Created data.tables")

# Add peptide length
colnames(PSP) <- c("protein", "peptide")
colnames(PCP) <- c("protein", "peptide")

Seq_stats$length <- Nmers
print("Added length")

### ---------------------------- (3) Save stats --------------------------------------
### Export protein stats
{
  Seq_stats_dir <- paste0(directory, "/Seq_stats/")
  suppressWarnings(dir.create(Seq_stats_dir))
}
vroom_write(Seq_stats, 
            delim = ",", num_threads = Ncpu,
            unlist(snakemake@output[["Seq_stats"]]))
print("Saved sequence stats")
print(Sys.time())


### ------------------------------------------ (4) Save all unique peptides ------------------------------------------
PCP %>%
  lazy_dt() %>%
  select(peptide) %>%
  mutate(index = str_sub(peptide, 1, index_length)) %>%
  group_by(index) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           pipe(sprintf("pigz > %s", paste0(directory, "/peptide_seqences/PCP_",.y$index, "_", filename, ".csv.gz"))),
                           delim = ",", num_threads = Ncpu))
print("Saved unique PCP")

PSP %>%
  lazy_dt() %>%
  select(peptide) %>%
  mutate(index = str_sub(peptide, 1, index_length)) %>%
  group_by(index) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           pipe(sprintf("pigz > %s", paste0(directory, "/peptide_seqences/PSP_",.y$index, "_", filename, ".csv.gz"))),
                           delim = ",", num_threads = Ncpu))
print("Saved unique PSP")


### ------------------------------------------ (5) Save Protein-peptide mapping ------------------------------------------
PCP %>%
  lazy_dt() %>%
  mutate(index = str_sub(peptide, 1, index_length)) %>%
  group_by(index) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           pipe(sprintf("pigz > %s", paste0(directory, "/peptide_mapping/PCP_map_",.y$index, "_", filename, ".csv.gz"))),
                           delim = ",", num_threads = Ncpu))
print("Saved PCP protein-peptide mapping")


PSP %>%
  lazy_dt() %>%
  mutate(index = str_sub(peptide, 1, index_length)) %>%
  group_by(index) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           pipe(sprintf("pigz > %s", paste0(directory, "/peptide_mapping/PSP_map_",.y$index, "_", filename, ".csv.gz"))),
                           delim = ",", num_threads = Ncpu))

print("Saved PSP protein-peptide mapping")
print(Sys.time())
