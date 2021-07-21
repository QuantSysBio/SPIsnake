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
dat = read.fasta(file="/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/inputSequence_1_2492.fasta", 
                 seqtype="AA", as.string = TRUE)
dat <- dat[1:100]

# Input proteome filter
min_protein_length = 9

# Nmers
Nmers = as.numeric(8:15)

# max intervening sequence length
MiSl = as.numeric(25)

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X)"

# Save into chunks according to first N letters
index_length = 1

### ---------------------------- (2) Compute PCP and PSP --------------------------------------
for (i in seq_along(Nmers)) {
  print(Sys.time())
  print(paste("Starting length: ",Nmers[i]))
  
  # Compute PCP and PSP per sequence
  dat_sort <- dat[names(sort(unlist(lapply(dat, nchar)), decreasing = T))]
  Pep_list <- mclapply(dat_sort, mc.cores = Ncpu, CutAndPaste_seq, nmer = Nmers[i], MiSl = MiSl)
  print(Sys.time())
  print("Computed PCP/PSP")
  
  # Tidy format
  PSP <- seq_list_to_dt(lapply(Pep_list, `[[`, 1))
  PCP <- seq_list_to_dt(lapply(Pep_list, `[[`, 2))
  Seq_stats <- rbindlist(lapply(Pep_list, `[[`, 3))
  print("Created data.tables")
  
  # Add peptide length
  colnames(PSP) <- c("protein", "peptide")
  colnames(PCP) <- c("protein", "peptide")
  # PCP$length <- Nmers[i]
  # PSP$length <- Nmers[i]
  Seq_stats$length <- Nmers[i]
  print("Added length")
  
  # Export protein-peptide mapping
  vroom_write(PCP, file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PCP_prot_pep_", Nmers[i], ".csv.gz"), 
              append = T, delim = ",", num_threads = Ncpu)
  print("Saved PCP")
  
  vroom_write(PSP, file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot_pep_", Nmers[i], ".csv.gz"), 
              append = T, delim = ",", num_threads = Ncpu)
  print("Saved PSP")
  
  # Export protein stats
  vroom_write(Seq_stats, file = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/Seq_stats.csv", 
              append = T, delim = ",", num_threads = Ncpu)
  print("Saved sequence stats")
  print(Sys.time())
  
  # Save all unique peptides
  PCP %>%
    lazy_dt() %>%
    select(peptide) %>%
    filter(!grepl(pattern = exclusion_pattern, peptide) | 
             grepl(pattern = "[:punct:]", peptide) |
             nchar(peptide) == 0) %>%
    mutate(index = substr(peptide, 1, index_length)) %>%
    group_by(index) %>%
    unique() %>%
    group_walk(~ vroom_write(.x, 
                             paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/unique_seqences/",.y$index, "_", Nmers[i], ".csv.gz"), 
                             append = T, delim = ",", num_threads = Ncpu))
  print("Saved unique PCP")
  
  PSP %>%
    lazy_dt() %>%
    select(peptide) %>%
    filter(!grepl(pattern = exclusion_pattern, peptide) | 
             grepl(pattern = "[:punct:]", peptide) |
             nchar(peptide) == 0) %>%
    mutate(index = substr(peptide, 1, index_length)) %>%
    group_by(index) %>%
    unique() %>%
    group_walk(~ vroom_write(.x, 
                             paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/unique_seqences/",.y$index, "_", Nmers[i], ".csv.gz"), 
                             append = T, delim = ",", num_threads = Ncpu))
  print("Saved unique PSP")
  rm(PCP, PSP, Seq_stats, Pep_list)
  gc()
  print(Sys.time())
}



# ----------------------------------------- DEV: ON -----------------------------------
a <- dat[names(sort(unlist(lapply(dat, nchar)), decreasing = T))]
# a <- mclapply(dat[1:50], mc.cores = Ncpu, CutAndPaste_seq, nmer = 8, MiSl = MiSl)


Seq_stats <- vroom("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/Seq_stats.csv", delim = ",",
                   col_names = c("protein","all_PCP","all_PSP","unique_PCP","unique_PSP","unique_PSP_noPCP","length"))


print(Sys.time())
Prot_pep <- vroom("/home/yhorokh/Snakemake/SPI-snake/long_tables/PSP_prot_pep.csv.gz", 
                  delim = ",",  altrep = T, num_threads = Ncpu, progress = T, 
                  col_names = c("protein", "peptide", "length"),
                  col_types = c("f", "c", "n")) %>%
  lazy_dt()
print(Sys.time())
pep_filt <- PCP %>%
  lazy_dt() %>%
  select(peptide) %>%
  left_join(Prot_pep)
print(Sys.time())


# Check saved file size difference
print("Start test")
for (j in 1:3) {
  vroom_write(PSP, file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot_pep_", Nmers[i], ".csv.gz"), 
              append = T, delim = ",", num_threads = Ncpu)
  print(paste("Saved PSP", j))
}


print("Start test 2")
PSP$protein <- stringr::str_split_fixed(PSP$protein, pattern = stringr::fixed("|"), n = 2)[,1]
for (j in 1:3) {
  vroom_write(PSP, file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot.split_pep_", Nmers[i], ".csv.gz"), 
              append = T, delim = ",", num_threads = Ncpu)
  print(paste("Saved PSP", j))
}




PSP_headers <- vroom(file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot_pep_8.csv.gz"), 
                     delim = ",",  altrep = T, num_threads = Ncpu, progress = T, 
                     col_names = c("protein", "peptide"),
                     col_types = c("f", "c"))
PSP_split <- vroom(file = paste0("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot.split_pep_8.csv.gz"), 
                     delim = ",",  altrep = T, num_threads = Ncpu, progress = T, 
                     col_names = c("protein", "peptide"),
                     col_types = c("f", "c"))



"/home/yhorokh/Snakemake/SPI-snake_backup_14.07.21/long_tables/PSP_prot_pep.csv.gz"






Join_pep_prot_chunk <- function(Peptide_dt, Prot_pep_map, chunk_size=10^6) {
  
  # Read protein-peptide mapping chunk-wise
  pep_out <- data.table(protein=NA, peptide=NA)
  start=1
  last_chunk=chunk_size
  chunk_size=chunk_size
  end=chunk_size
  while (last_chunk == chunk_size) {
    print(paste("Processing rows:", start, "-", end))
    
    
    ### Reading chunked file is impossible with vroom
    Prot_pep <- vroom(
      file = Prot_pep_map,
      skip = start - 1,
      n_max = chunk_size,
      delim = ",",
      altrep = T,
      col_names = c("protein", "peptide"),
      col_types = c("f", "c"),
      # progress = T,
      num_threads = Ncpu
      ) %>%
      as.data.table()
    
    
    ### Reading compressed file is slow with fread
    # Prot_pep <- fread(file = Prot_pep_map, 
    #                   sep = ",", 
    #                   skip = start - 1,
    #                   nrows = chunk_size, 
    #                   nThread=Ncpu,
    #                   header = F,
    #                   # col.names = c("protein", "peptide"),
    #                   # select=c(protein="factor", peptide="character"),
    #                   showProgress=TRUE,
    #                   data.table=TRUE) %>%
    #   lazy_dt() %>%
    #   rename(protein=V1, peptide=V2)
    
    pep_filt <- Peptide_dt %>%
      lazy_dt() %>%
      select(peptide) %>%
      left_join(Prot_pep) %>%
      filter(!is.na(protein)) %>% 
      as.data.table()
    pep_out <- rbind(pep_out, pep_filt)
    
    # Update next chunk read-in
    last_chunk <- nrow(as.data.table(Prot_pep))
    end <- end + last_chunk
    start <- start + last_chunk
  }
  pep_out <- pep_out %>%
    filter(!is.na(protein)) 
  return(pep_out)
}

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1)
# Sys.setenv("VROOM_CONNECTION_SIZE" = 10^5)
test <- Join_pep_prot_chunk(Peptide_dt = PCP, 
                            Prot_pep_map = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot_pep_8.csv.gz", 
                            chunk_size = 10^3)




### Arrow
library(arrow)
Prot_pep <- read_csv_arrow(Prot_pep_map, 
                           col_names = c("protein", "peptide"),
                           as_data_frame=F,
)




pep_filt <- Peptide_dt %>%
  lazy_dt() %>%
  select(peptide) %>%
  left_join(Prot_pep) %>%
  filter(!is.na(protein)) %>% 
  as.data.table()
pep_out <- rbind(pep_out, pep_filt)

# Compare parquet vs csv for writing
print(Sys.time())
write_parquet(PSP, 
              sink = "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/parquet", 
              compression = "snappy")
print(Sys.time())
vroom_write(PSP, "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/PSP_prot_pep_8_test.csv.gz", 
            append = T, delim = ",", num_threads = Ncpu)
print(Sys.time())



print(Sys.time())
write_dataset(PSP, "/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/parquet", format = "parquet")
print(Sys.time())


print(Sys.time())

a <- Table$create(PSP)

print(Sys.time())
PSP %>%
  lazy_dt() %>%
  select(peptide) %>%
  filter(!grepl(pattern = exclusion_pattern, peptide) | 
           grepl(pattern = "[:punct:]", peptide) |
           nchar(peptide) == 0) %>%
  mutate(index = substr(peptide, 1, index_length)) %>%
  Table$create() %>%
  group_by(index) %>%
  write_dataset("/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/parquet", 
                format = "parquet", 
                compression = "snappy")
print(Sys.time())

