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

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

# Manual startup
# {
#   ### setwd("/home/yhorokh/Snakemake/SPIsnake-main")
#   directory = "results/DB_exhaustive/"
#   dir_DB_Fasta_chunks = "results/DB_exhaustive/Fasta_chunks/"
#   Master_table_expanded <- read.csv("results/DB_exhaustive/Master_table_expanded.csv")
#   filename = "results/DB_exhaustive/Seq_stats/15_cis-PSP_Minimal_25_SwissProt_UP000005640_11600_12273.fasta.csv.gz"
#   index_length = 2
#   max_protein_length = 100
#   Ncpu = 7
# 
#   filename <- str_remove(filename, ".csv.gz") %>%
#     str_split_fixed(pattern = fixed("Seq_stats/"), n = 2)
#   filename <- filename[,2]
#   print(filename)
# 
#   params <- Master_table_expanded[Master_table_expanded$filename == filename,]
#   print(t(params))
# }

### ---------------------------- (1) Read input file and extract info ----------------------------
# Master_table_expanded
Master_table_expanded <- read.csv(snakemake@input[["Master_table_expanded"]])

# Fasta input
dir_DB_Fasta_chunks = snakemake@params[["dir_DB_Fasta_chunks"]]

# Output dir
directory = snakemake@params[["directory"]]
suppressWarnings(dir.create(paste0(directory, "/peptide_seqences")))
suppressWarnings(dir.create(paste0(directory, "/peptide_mapping")))

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
# Ncpu = availableCores(methods = "Slurm")
Ncpu = snakemake@params[["cpus_for_R"]]
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl, debug=T)
data.table::setDTthreads(Ncpu)

print(paste0("number of CPUs: ", Ncpu))

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

# fst compression level
fst_compression = as.integer(snakemake@params[["fst_compression"]])

### ---------------------------- (2) Compute PCP and PSP --------------------------------------
print(Sys.time())
print(paste("Starting length: ", Nmers))

# Compute PCP and PSP per sequence
dat_sort <- dat[names(sort(unlist(lapply(dat, nchar)), decreasing = T))]

if (grepl("cis-PSP", Splice_type)==T) {
  print("Computing PCP and cis-PSP")
  
  Pep_list <- bettermc::mclapply(dat_sort, 
                       mc.cores = Ncpu,
                       mc.cleanup=T, mc.preschedule=F, mc.retry = 3,
                       CutAndPaste_seq_from_big_sp_fast, 
                       big_sp_input=index_list_result,
                       nmer = Nmers, 
                       MiSl=MiSl)
  print(Sys.time())
  print("Computed PCP/PSP")
} else if (Splice_type == "PCP") {
  print("Computing PCP only")
  
  Pep_list <- bettermc::mclapply(dat_sort, mc.cores = Ncpu,  mc.cleanup=T, mc.preschedule=F,
                        CutAndPaste_seq_PCP, nmer = Nmers)
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

### ------------------------------------------ (4) Save peptides and mapping ------------------------------------------
PCP[, index := str_sub(peptide, start = 1, end = index_length)] %>%
  split(by = "index", drop = T, keep.by = T) %>%
  bettermc::mclapply(mc.cores = Ncpu, mc.cleanup=T, mc.preschedule=F, FUN = function(x){
    x <- x[!(str_detect(peptide, exclusion_pattern) | !str_length(peptide) == Nmers)]
    
    if (nrow(x) > 0) {
      index <- x$index[[1]]
      
      x[, .(protein, peptide)] %>% 
        write_fst(path = paste0(directory, "/peptide_mapping/PCP_map_", index, "_", filename, ".fst"), compress = fst_compression) 
      
      x[, .(peptide)] %>% 
        unique() %>%
        write_fst(path = paste0(directory, "/peptide_seqences/PCP_", index, "_", filename, ".fst"), compress = fst_compression) 
    }
  }) %>%
  capture.output()
print("Saved PCP")

PSP[, index := str_sub(peptide, start = 1, end = index_length)] %>%
  split(by = "index", drop = T, keep.by = T) %>%
  bettermc::mclapply(mc.cores = Ncpu, mc.cleanup=T, mc.preschedule=F, FUN = function(x){
    x <- x[!(str_detect(peptide, exclusion_pattern) | !str_length(peptide) == Nmers)]
    
    if (nrow(x) > 0) {
      index <- x$index[[1]]
      
      x[, .(protein, peptide)] %>% 
        write_fst(path = paste0(directory, "/peptide_mapping/PSP_map_", index, "_", filename, ".fst"), compress = fst_compression) 
      
      x[, .(peptide)] %>% 
        unique() %>%
        write_fst(path = paste0(directory, "/peptide_seqences/PSP_", index, "_", filename, ".fst"), compress = fst_compression)
    }
  }) %>%
  capture.output()
print("Saved PSP")
print(Sys.time())

print("----- memory usage by Slurm -----")
jobid = system("echo $SLURM_JOB_ID")
system(paste0("sstat ", jobid)) %>%
  print()

system("sacct --format='JobID,JobName,State,Elapsed,AllocNodes,NCPUS,NodeList,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask,ReqMem,MaxDiskWrite'") %>%
  print()

print("----- memory usage by R -----")
memory.profile() %>%
  print()

print("----- connections -----")
showConnections() %>%
  print()

print("----- removing cluster -----")
print(cl)
parallel::stopCluster(cl)

print("----- garbage collection -----")
gc() %>%
  print()
