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
sink(log, split = TRUE)

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cleaver))
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
seq_vectorized <- Vectorize(seq.default, vectorize.args = c("from", "to"), SIMPLIFY = F) 
print(sessionInfo())

### Create temporary directory for vroom
{
  Sys.getenv("TMPDIR") %>% print()
  Sys.getenv("VROOM_TEMP_PATH") %>% print()
  
  vroom_dir = "/tmp/vroom"
  suppressWarnings(dir.create(vroom_dir))
  Sys.setenv(VROOM_TEMP_PATH = vroom_dir)
  Sys.getenv("VROOM_TEMP_PATH") %>%print()
  
  tmp_file = tempfile()
  print(tmp_file)
}

# {
#   ### Manual startup
#   ### setwd("/home/yhorokh/Snakemake/SPIsnake-main")
#   directory = "results/DB_exhaustive/"
#   dir_DB_Fasta_chunks = "results/DB_exhaustive/Fasta_chunks/"
#   suppressWarnings(dir.create(paste0(directory, "/peptide_seqences")))
#   suppressWarnings(dir.create(paste0(directory, "/peptide_mapping")))
#   Master_table_expanded <- read.csv("results/DB_exhaustive/Master_table_expanded.csv")
#   filename = "results/DB_exhaustive/Seq_stats/5_30_trypsin_2_25_lncRNA_2395202_4790402_1.fasta.csv.gz"
#   index_length = 1
#   max_protein_length = 500
#   exclusion_pattern <- "(U|X|O|\\*)"
#   fst_compression = 100
# replace_I_with_L = TRUE
#   
#   Ncpu = 63
#   cl <- parallel::makeForkCluster(Ncpu)
#   cl <- parallelly::autoStopCluster(cl, debug=T)
#   data.table::setDTthreads(Ncpu)
#   
#   filename <- str_remove(filename, ".csv.gz") %>%
#     str_split_fixed(pattern = fixed("Seq_stats/"), n = 2)
#   filename <- filename[,2]
#   print(filename)
#   
#   params <- Master_table_expanded[Master_table_expanded$filename == filename,] %>%
#     select(Proteome, Splice_type, N_mers, Min_Interv_length, chunk) %>%
#     unique()
#   print(t(params))
#   Splice_type = as.character(params$Splice_type)
#   MiSl = as.numeric(params$Min_Interv_length)
#   # Nmers
#   if (Splice_type == "PCP" | (grepl("cis-PSP", Splice_type) == TRUE)) {
#     # For PCP and cis-PSP corresponds to Nmer to be produced
#     Nmers = as.numeric(params$N_mers)
#   } else {
#     # For cleaver digestions defines min and max peptide length
#     Nmers = as.numeric(unlist(str_split(params$N_mers, "_", 2)))
#   }
#   
#   index_list_result = paste(Nmers, MiSl, sep = "_")
#   if (grepl("cis-PSP", Splice_type)==T) {
#     load(paste0(directory, "/PSP_indices/", index_list_result, ".rds"))
#   }
#   
#   proteome = list.files(dir_DB_Fasta_chunks, pattern = params$chunk, recursive = T, full.names = T)
#   proteome = proteome[str_ends(proteome, ".fasta")]
#   dat = readAAStringSet(proteome)
#   
#   # Keep only proteome name
#   proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
#   proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]
#   enzymes <- c("arg-c proteinase", "asp-n endopeptidase", "bnps-skatole", "caspase1", "caspase2", "caspase3", "caspase4", "caspase5", "caspase6", "caspase7", "caspase8", "caspase9", "caspase10", "chymotrypsin-high", "chymotrypsin-low", "clostripain", "cnbr", "enterokinase", "factor xa", "formic acid", "glutamyl endopeptidase", "granzyme-b", "hydroxylamine", "iodosobenzoic acid", "lysc", "lysn", "neutrophil elastase", "ntcb", "pepsin1.3", "pepsin", "proline endopeptidase", "proteinase k", "staphylococcal peptidase i", "thermolysin", "thrombin", "trypsin")
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

params <- Master_table_expanded[Master_table_expanded$filename == filename,] %>%
  select(Proteome, Splice_type, N_mers, Min_Interv_length, chunk) %>%
  unique()
print(t(params))

### Extract parameters
# Input fasta
proteome = list.files(dir_DB_Fasta_chunks, pattern = params$chunk, recursive = T, full.names = T)
proteome = proteome[str_ends(proteome, ".fasta")]
dat = readAAStringSet(proteome)
# if (replace_I_with_L == TRUE) {
#   dat <- chartr("I", "L", dat)
# }

# Keep only proteome name
proteome <- unlist(strsplit(proteome, "/", fixed = T))[grep(".fasta", unlist(strsplit(proteome, "/", fixed = T)))]
proteome <- unlist(strsplit(proteome, ".fasta", fixed = T))[1]

# Splice type
Splice_type = as.character(params$Splice_type)

# max intervening sequence length
MiSl = as.numeric(params$Min_Interv_length)

# Nmers
if (Splice_type == "PCP" | (grepl("cis-PSP", Splice_type) == TRUE)) {
  # For PCP and cis-PSP corresponds to Nmer to be produced
  Nmers = as.numeric(params$N_mers)
} else {
  # For cleaver digestions defines min and max peptide length
  Nmers = as.numeric(unlist(str_split(params$N_mers, "_", 2)))
}

# CPUs
Ncpu = snakemake@params[["cpus_for_R"]]
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl, debug=T)
data.table::setDTthreads(Ncpu)

print(paste0("number of CPUs: ", Ncpu))

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X|O|\\*)"

# Save into chunks according to first N letters
index_length = as.integer(snakemake@params[["AA_index_length"]])

# Save into chunks according to first N letters
max_protein_length = as.integer(snakemake@params[["max_protein_length"]])

# # Replace I with L
# replace_I_with_L = as.logical(snakemake@params[["replace_I_with_L"]])
# print(paste("Replace I with L:", replace_I_with_L))

# Load pre-computed index to speed-up PSP generation
index_list_result = paste(Nmers, MiSl, sep = "_")
if (grepl("cis-PSP", Splice_type)==T) {
  load(paste0(directory, "/PSP_indices/", index_list_result, ".rds"))
}

# fst compression level
fst_compression = as.integer(snakemake@params[["fst_compression"]])

# Cleaver params
enzymes <- c("arg-c proteinase", "asp-n endopeptidase", "bnps-skatole", "caspase1", "caspase2", "caspase3", "caspase4", "caspase5", "caspase6", "caspase7", "caspase8", "caspase9", "caspase10", "chymotrypsin-high", "chymotrypsin-low", "clostripain", "cnbr", "enterokinase", "factor xa", "formic acid", "glutamyl endopeptidase", "granzyme-b", "hydroxylamine", "iodosobenzoic acid", "lysc", "lysn", "neutrophil elastase", "ntcb", "pepsin1.3", "pepsin", "proline endopeptidase", "proteinase k", "staphylococcal peptidase i", "thermolysin", "thrombin", "trypsin")
custom_trypsin = c("[KR](?=\\w)")

### ---------------------------- (2) Cleave PCP --------------------------------------
print(Sys.time())
dat_sort <- dat[width(dat) >= min(Nmers)]

if (Splice_type == "PCP" | (grepl("cis-PSP", Splice_type) == TRUE)) {
  
  # Pre-filter and sort input sequences
  print(paste("Starting length: ", Nmers))
  dat_sort <- dat_sort[order(width(dat_sort), decreasing = T),]
  print("Computing PCP")
  
  split_chunks <- rep(1:(5*Ncpu), length.out=length(dat_sort))
  PCP <- split(dat_sort, split_chunks) %>%
    bettermc::mclapply(mc.preschedule = TRUE, 
                       mc.cores = Ncpu, 
                       mc.cleanup = TRUE, 
                       mc.force.fork = TRUE, 
                       mc.retry = 3, 
                       Nmers = Nmers,
                       FUN = function(x, Nmers){
                         extractAt(x, IRangesList(start = seq_vectorized(from = 1, to = width(x)-Nmers+1),
                                                  end = seq_vectorized(from = Nmers, to = width(x))))
                       })
  # Tidy format
  PCP <- rbindlist(lapply(PCP, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  setkey(PCP, peptide, protein)
  PCP <- PCP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein), by=index] %>%
    .[, type := "PCP"]
  
  # gc()
  print(Sys.time())
  print("Computed PCP")
  
  Seq_stats <- PCP
} 

### ---------------------------- (3) Splice PSP --------------------------------------
if (grepl("cis-PSP", Splice_type) == TRUE) {
  print("Computing cis-PSP")
  
  # Match indices to proteins of corresponding length
  dat_sort <- dat_sort[width(dat_sort) > Nmers]
  inputs <- vector(mode = "list", length = length(dat_sort))
  for (i in seq_along(inputs)) {
    inputs[[i]][[1]] <- dat_sort[[i]]
    inputs[[i]][[2]] <- index_list_result[[width(dat_sort[i])]]
  }
  
  # Generate PSP sequences
  PSP <- bettermc::mclapply(X = inputs, 
                            Nmers = Nmers,
                            FUN = Generate_PSP_2,
                            mc.cores = Ncpu,
                            mc.cleanup = TRUE, 
                            mc.preschedule = TRUE, 
                            mc.force.fork = TRUE, 
                            mc.retry = 3)
  names(PSP) <- names(dat_sort)
  PSP <- setDT(stack(PSP))
  setnames(PSP, c("peptide", "protein"))
  setkey(PSP, peptide, protein)
  PSP <- PSP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein), by=index] %>%
    .[, type := "PSP"]
  
  print(Sys.time())
  print("Computed PCP/PSP")
  
  Seq_stats <- rbindlist(list(PCP, PSP))
} 

### ---------------------------- (4) Cleave enzymatic digestions --------------------------------------
if (str_detect(Splice_type, str_c(enzymes, collapse = "|"))) {
  
  enzym <- str_split_fixed(Splice_type, "_", 2)[,1]
  missedCleavages <- str_split_fixed(Splice_type, "_", 2)[,2]
  
  print(paste("Computing digestions by", enzym, "with", missedCleavages, "missed cleavages"))
  print(paste("Peptides of length:", min(Nmers), "-", max(Nmers), "will be kept"))
  
  # Generate in silico enzymatic digestions
  dat_sort <- dat[width(dat) >= min(Nmers)]
  
  ### Add the same sequences without N-terminal Met
  {
    add_X <- dat_sort[str_starts(dat_sort, "M|X")]
    add_X <- subseq(add_X, start = 2, end = ifelse(width(add_X) >= 2*max(Nmers), 2*max(Nmers), width(add_X)))
    dat_sort <- c(dat_sort, add_X)
  }
  
  split_chunks <- rep(1:(5*Ncpu), length.out=length(dat_sort))
  if (enzym == "trypsin") {
    enzym_out <- split(dat_sort, split_chunks) %>%
      bettermc::mclapply(FUN = cleave,
                         missedCleavages = 0:missedCleavages,
                         custom = custom_trypsin, 
                         unique = FALSE,
                         mc.cores = Ncpu,
                         mc.cleanup = TRUE, 
                         mc.preschedule = TRUE, 
                         mc.force.fork = TRUE, 
                         mc.retry = 3)
  } else {
    enzym_out <- split(dat_sort, split_chunks) %>%
      bettermc::mclapply(FUN = cleave,
                         enzym = enzym, 
                         missedCleavages = 0:missedCleavages,
                         unique = FALSE,
                         mc.cores = Ncpu,
                         mc.cleanup = TRUE, 
                         mc.preschedule = TRUE, 
                         mc.force.fork = TRUE, 
                         mc.retry = 3)
  }
  
  # Tidy format
  enzym_out <- rbindlist(lapply(enzym_out, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  setkey(enzym_out, peptide, protein)
  enzym_out <- enzym_out %>% 
    unique() %>%
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein, index), by=index] %>%
    .[,length := str_length(peptide), by=index] %>%
    .[inrange(length, lower = min(Nmers), upper = max(Nmers), incbounds=TRUE), .(peptide, length, protein), by=index] %>%
    .[, type := Splice_type] 
  
  print(Sys.time())
  print(paste("Computed", Splice_type))
  
  Seq_stats <- enzym_out
  if (nrow(enzym_out) == 0) {
    rm(enzym_out)
  }
}

### ---------------------------- (5) Save stats --------------------------------------
if (Splice_type == "PCP" | (grepl("cis-PSP", Splice_type) == TRUE)) {
  # Add peptide length
  setkey(Seq_stats, protein, type)
  Seq_stats <- Seq_stats %>%
    lazy_dt() %>%
    group_by(protein, type) %>%
    summarise(total_peptides = n(),
              unique_peptides = n_distinct(peptide)) %>%
    mutate(length = Nmers) %>%
    as.data.table()
  
  print("Added length")
} else if (str_detect(Splice_type, str_c(enzymes, collapse = "|"))) {
  # Add peptide length
  setkey(Seq_stats, protein, type, length)
  Seq_stats <- Seq_stats %>%
    lazy_dt() %>%
    group_by(protein, type, length) %>%
    summarise(total_peptides = n(),
              unique_peptides = n_distinct(peptide)) %>%
    as.data.table()
  
  print("Added length")
}

Seq_stats_dir <- paste0(directory, "/Seq_stats/")
suppressWarnings(dir.create(Seq_stats_dir))

fwrite(Seq_stats,
       sep = ",", nThread = Ncpu,
       file = unlist(snakemake@output[["Seq_stats"]]))
print("Saved sequence stats")
print(Sys.time())

### ------------------------------------------ (6) Save peptides and mapping ------------------------------------------
data.table::setDTthreads(Ncpu)
if (exists("PCP")) {
  PCP %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_mapping/PCP_map_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("protein", "peptide")] %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_seqences/PCP_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("peptide")]
  print("Saved PCP")
  print(Sys.time())
}
if (exists("PSP")) {
  PSP %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_mapping/PSP_map_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("protein", "peptide")] %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_seqences/PSP_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("peptide")]
  print("Saved PSP")
  print(Sys.time())
}
if (exists("enzym_out")) {
  enzym_out %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_mapping/", Splice_type, "_map_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("protein", "peptide")] %>% 
    .[, write_fst(unique(.SD), paste0(directory, "/peptide_seqences/", Splice_type, "_", .BY, "_", filename, ".fst"), compress = fst_compression), 
      by=index, .SDcols=c("peptide")]
  print("Saved enzym_out")
  print(Sys.time())
}

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
sink()
