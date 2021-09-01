library(data.table)
library(dtplyr)
library(dplyr)
library(seqinr)
library(parallel)
library(parallelly)
library(foreach)
library(vroom)
source("/home/yhorokh/Snakemake/SPI-snake/src/snakefiles/functions.R")

# Functions
seq_list_to_dt <- function(seq_list){
  rbindlist(lapply(seq_list, as.data.table), idcol = "id")
}

### ---------------------------- (1) Read input file and extract info ----------------------------
# Input fasta
proteome = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/fiveUTR_seq_3_frame/fiveUTR_seq_3_frame_1_2800.fasta"
dat = read.fasta(file=proteome, 
                 seqtype="AA", as.string = TRUE)

# Keep only proteome name
proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]

# Input proteome filter
min_protein_length = 9

# Nmers
Nmers = as.numeric(15)

# max intervening sequence length
MiSl = as.numeric(25)

# Splice type
Splice_type = "PCP, cis-PSP"

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X|\\*)"

# Save into chunks according to first N letters
index_length = 1

# Temp - Snakemake will create directories automatically
dir.create("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_seqences/")
dir.create("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_mapping/")

# Output dir
directory = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive"
suppressWarnings(
 dir.create(paste0(directory, "/peptide_seqences"))
)
suppressWarnings(
  dir.create(paste0(directory, "/peptide_mapping"))
)

### ---------------------------- (2) Compute PCP and PSP --------------------------------------
print(Sys.time())
print(paste("Starting length: ",Nmers))

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

### ---------------------------- (3) Export --------------------------------------
### Export protein stats
vroom_write(Seq_stats, file = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/Seq_stats.csv", 
            append = T, delim = ",", num_threads = Ncpu)
print("Saved sequence stats")
print(Sys.time())

### Save all unique peptides
PCP %>%
  lazy_dt() %>%
  select(peptide) %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
          nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_seqences/",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
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
                           paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_seqences/",.y$index, "_", Nmers, "_", proteome, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved unique PSP")


### Save Protein-peptide mapping
PCP %>%
  lazy_dt() %>%
  filter(!(grepl(pattern = exclusion_pattern, peptide) | 
             nchar(peptide) == 0)) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  group_by(index) %>%
  unique() %>%
  group_walk(~ vroom_write(.x, 
                           paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_mapping/PCP_map_", .y$index, "_", Nmers, ".csv.gz"), 
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
                           paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/peptide_mapping/PSP_map_", .y$index, "_", Nmers, ".csv.gz"), 
                           append = T, delim = ",", num_threads = Ncpu))
print("Saved PSP protein-peptide mapping")
print(Sys.time())


