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
# author:       YH, JL

library(data.table)
library(dtplyr)
library(dplyr)
library(seqinr)
library(parallel)
library(parallelly)
library(foreach)
library(vroom)
source("src/snakefiles/functions.R")

# Functions
seq_list_to_dt <- function(seq_list){
  rbindlist(lapply(seq_list, as.data.table), idcol = "id")
}

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master_table_expanded
Master_table <- read.csv(snakemake@input[["Master_table_expanded"]])
print(Master_table)

# Output dir
directory = snakemake@params[["directory"]]
# directory = "/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive"
suppressWarnings(
  dir.create(paste0(directory, "/peptide_seqences"))
)
suppressWarnings(
  dir.create(paste0(directory, "/peptide_mapping"))
)

# Wildcard
# filename = "results/DB_exhaustive/.Generate_cis-PSP_15_25_Unimod_proteome_expressed_gencode_2177_2722.fasta.done"
filename = snakemake@output[[1]]
filename <- str_remove(filename, ".done") %>%
  str_split_fixed(pattern = fixed("/."), n = 2)
filename <- filename[,2] %>%
  str_remove(pattern = "Generate_")
print(filename)

params <- Master_table_expanded[Master_table_expanded$filename == filename,]
print(t(params))

### Extract parameters
# Input fasta
proteome = list.files(directory, pattern = params$chunk, recursive = T, full.names = T)
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

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X|\\*)"

# Save into chunks according to first N letters
index_length = 1

### ---------------------------- (2) Compute PCP and PSP --------------------------------------
print(Sys.time())
print(paste("Starting length: ", Nmers))

# Compute PCP and PSP per sequence
dat_sort <- dat[names(sort(unlist(lapply(dat, nchar)), decreasing = T))]

if (grepl("cis-PSP", Splice_type)==T) {
  print("Computing PCP and cis-PSP")
  
  Pep_list <- mclapply(dat_sort, mc.cores = Ncpu, CutAndPaste_seq, nmer = Nmers, MiSl = MiSl)
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
            pipe(sprintf("pigz > %s", paste0(Seq_stats_dir, filename, ".csv.gz"))))
print("Saved sequence stats")
print(Sys.time())


### ------------------------------------------ (4) Save all unique peptides ------------------------------------------
PCP %>%
  lazy_dt() %>%
  select(peptide) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive/peptide_seqences/PCP_",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved unique PCP")

PSP %>%
  lazy_dt() %>%
  select(peptide) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive/peptide_seqences/PSP_",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved unique PSP")

### ------------------------------------------ (5) Save Protein-peptide mapping ------------------------------------------
PCP %>%
  lazy_dt() %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive/peptide_mapping/PCP_map_",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved PCP protein-peptide mapping")

PSP %>%
  lazy_dt() %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake_dev/results/DB_exhaustive/peptide_mapping/PSP_map_",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved PSP protein-peptide mapping")
print(Sys.time())

### ------------------------------------------- Bookmark ------------------------------------------------------------------------