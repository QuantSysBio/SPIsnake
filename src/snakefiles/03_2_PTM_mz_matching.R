### ---------------------------------------------- Aggregate_PTM_mz_matching ----------------------------------------------
# description:  Find a unique set of peptide sequences, compute molecular weight (MW) and do m/z matching with all input mass_lists. 
#               If PTMs are required, generate them too. 
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .csv.gz
#               - Unique peptides MW .csv.gz
#               - Unique peptides after m/z matching per experiment .csv.gz
#               - Peptide-mass_list matching
#               - Peptide filtering stats .csv
#               
# author:       YH, JL, KP

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(data.table))

source("src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
{
  # # Manual setup
  # Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv")
  # Peptide_aggregation_table <- vroom("results/DB_PTM_mz/Peptide_aggregation_table.csv", delim = ",")
  # Experiment_design <- vroom("data/Experiment_design.csv", delim = ",")
  # dir_DB_exhaustive = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_exhaustive"
  # dir_DB_PTM_mz = "/home/yhorokh/SNAKEMAKE/SPIsnake-main/results/DB_PTM_mz"
  # suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory")))
  # 
  # filename = "results/DB_PTM_mz/chunk_aggregation_status/F_9.csv"
  # filename <- filename %>%
  #  str_split_fixed(pattern = fixed("chunk_aggregation_status/"), n = 2)
  # filename <- filename[,2] %>%
  #  str_remove(pattern = ".csv")
  # print(filename)
  # MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
  #  as_tibble() %>%
  #  mutate(file = str_remove_all(value, ".txt"))
  # 
  # ### CPUs
  # Ncpu = availableCores()
  # cl <- parallel::makeForkCluster(Ncpu)
  # cl <- parallelly::autoStopCluster(cl)
  # setDTthreads(Ncpu)
  # 
  # # Save into chunks according to first N letters
  # index_length = 1
}

# Experiment_design
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)

# Master_table_expanded
Master_table_expanded <- vroom(snakemake@input[["Master_table_expanded"]], show_col_types = FALSE)

# Master_table_expanded
Peptide_aggregation_table <- vroom(snakemake@input[["Peptide_aggregation_table"]], delim = ",", show_col_types = FALSE)

# Output dir
dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory")))

### Wildcard
filename = snakemake@output[[1]]
filename <- filename %>%
  str_split_fixed(pattern = fixed("chunk_aggregation_status/"), n = 2)
filename <- filename[,2] %>%
  str_remove(pattern = ".csv")
print(filename)

# Save into chunks according to first N letters
index_length = as.integer(snakemake@params[["AA_index_length"]])

### Mass_lists
# Only user defined lists will be used 
MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
  as_tibble() %>%
  mutate(mass_list = str_remove_all(value, ".txt")) %>%
  filter(mass_list %in% Experiment_design$Filename) %>%
  rename(mass_list_file = value) %>%
  mutate(AA_length = filename) 

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (2) Operation mode --------------------------------------
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/files_mz_match/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/")))

# Check if there exist previous outputs to be updated:
processed_files <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory"), pattern = paste0(filename, ".csv"))

operation_mode = ifelse(length(processed_files) > 0, "Update", "Generation")
print(paste0("Mode: ", operation_mode))

if (operation_mode == "Update") {
  print("Updating file list to be processed")
  processed_files <- vroom(file = paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/", processed_files), delim = ",",  show_col_types = FALSE)
  # processed_files <- processed_files[1:60,]
}

### ---------------------------- (3) Uniqueness --------------------------------------
# Select all files with the same {AA-index}{length}
peptide_chunks <- list.files(paste0(dir_DB_exhaustive, "/peptide_seqences"), pattern = ".csv.gz", full.names = T) %>%
  as_tibble() %>%
  mutate(file = str_remove_all(value, ".csv.gz")) %>%
  mutate(file = str_split_fixed(file, "/peptide_seqences/", n = 2)[,2]) %>%
  mutate(Splice_type = ifelse(str_starts(file, "PCP_"), "PCP", "other")) %>%
  mutate(Splice_type = ifelse(str_starts(file, "PSP_"), "PSP", Splice_type)) %>%
  mutate(file = str_sub(file, start=5)) %>%
  mutate(file = str_remove(file, "PCP_")) %>%
  mutate(file = str_remove(file, "cis-PSP_")) %>%
  filter(str_starts(file, filename)) %>%
  mutate(AA_length = filename) 

if (operation_mode == "Update") {
  # Don't use absolute path when checking for file completeness
  tmp1 <- peptide_chunks %>%
    select(-value) %>%
    left_join(MS_mass_lists)
  tmp2 <- select(processed_files, colnames(tmp1))
  keep <- anti_join(tmp1, tmp2)
  
  # For new peptide sequences
  peptide_chunks <- peptide_chunks %>%
    filter(file %in% keep$file & 
             Splice_type %in% keep$Splice_type)
  
  # For new mass lists
  MS_mass_lists <- MS_mass_lists %>%
    filter(mass_list %in% keep$mass_list)
  
  rm(tmp1, tmp2, keep)
}

# End the script if no files to be updated
if (nrow(peptide_chunks) == 0) {
  print("Looks like there are no new peptides to be processed")
  
  save_peptide_chunks <- peptide_chunks %>%
    rbind(processed_files) %>%
    mutate(Time = Sys.time()) 
  
  vroom_write(save_peptide_chunks, delim = ",", num_threads = Ncpu,
              unlist(snakemake@output[["chunk_aggregation_status"]]))
} else {
  
  peptides <- as.list(peptide_chunks$value)
  names(peptides) <- peptide_chunks$file
  head(peptides)
  # peptides <- peptides[1:10]
  
  # PCP
  if (length(peptides[str_detect(peptides, "/PCP_")]) > 0) {
    PCP_list <- peptides[str_detect(peptides, "/PCP_")]  %>%
      lapply(FUN = vroom, delim = ",", num_threads = Ncpu, show_col_types = FALSE) 
  }
  # PSP
  if (length(peptides[str_detect(peptides, "/PSP_")]) > 0) {
    PSP_list <- peptides[str_detect(peptides, "/PSP_")]  %>%
      mclapply(FUN = vroom, delim = ",", num_threads = Ncpu, show_col_types = FALSE)
  }
  
  ### ---------------------------- (4) Compute MW --------------------------------------
  ### And check for peptide uniqueness
  pep_types <- c()
  if (length(PCP_list) > 0) {
    PCP <- PCP_list %>%
      rbindlist() %>%
      lazy_dt() %>%
      unique() %>%
      mutate(MW=computeMZ_biostrings(peptide)) %>%
      as.data.table()
    
    pep_types <- c(pep_types, "PCP")
  }
  
  if (length(PSP_list) > 0) {
    PSP <- PSP_list %>%
      rbindlist() %>%
      lazy_dt() %>%
      select(peptide) %>%
      unique() %>%
      mutate(MW=computeMZ_biostrings(peptide)) %>%
      as.data.table()
    
    pep_types <- c(pep_types, "PSP")
  }
  
  ### ---------------------------- (5) m/z matching --------------------------------------
  mz_nomod <- list()
  mz_stats <- tibble()
  for (i in pep_types) {
    input <- get(i)
    
    for (j in 1:nrow(MS_mass_lists)) {
      print(MS_mass_lists$mass_list[j])
      MS_mass_list <- MS_mass_lists$mass_list[j]
      
      tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
      mzList = vroom(file = paste0("data/MS_mass_lists/", MS_mass_lists$mass_list_file[j]), 
                     delim = ",", num_threads = Ncpu, 
                     col_names = "Precursor_mass", show_col_types = FALSE) %>%
        lazy_dt() %>%
        mutate(Min = Precursor_mass - Precursor_mass * tolerance * 10 ** (-6)) %>%
        mutate(Max = Precursor_mass + Precursor_mass * tolerance * 10 ** (-6)) %>%
        select(-Precursor_mass)  %>%
        unique() %>%
        as.data.table()
      
      ### Which peptide sequences pass the MW filter
      mz_nomod[[i]][[MS_mass_list]] <- input[mixANDmatch3(mzMin=mzList$Min, mzMax=mzList$Max, MW0 = input$MW),] %>%
        na.omit()
      
      ### Save stats
      mz_stats_ij = tibble(Splice_type=i,
                           mass_list = MS_mass_list,
                           mz_matched_peptides = nrow(mz_nomod[[i]][[MS_mass_list]]))
      mz_stats = rbind(mz_stats, mz_stats_ij)
    }
  }
  
  ### ---------------------------- (6) Generate PTMs --------------------------------------
  # to be implemented by @Kamil
  
  
  
  
  ### ---------------------------- (7) Peptide sequence export --------------------------------------
  # Unmodified m/z matched sequences
  for (i in seq_along(mz_nomod)) {
    for (j in seq_along(mz_nomod[[i]])) {
      mz_nomod[[i]][[j]] %>%
        vroom_write(paste0(dir_DB_PTM_mz, "/files_mz_match/", names(mz_nomod)[i], "_", filename, "_", names(mz_nomod[[i]])[j], ".csv.gz"),
                    delim = ",", num_threads = Ncpu, append = TRUE)
    }
  }
  
  ### --------------------------------------(8) Chunk aggregation status (Snakemake output) --------------------------------------
  ### Record aggregation stats
  {
    aggregation_stats <- tibble(total_PCP = sum(unlist(lapply(PCP_list, nrow))), 
                                unique_PCP = nrow(PCP),
                                total_PSP = sum(unlist(lapply(PSP_list, nrow))),
                                unique_PSP = nrow(PSP)
    )
    print("Loaded peptides")
    print(t(aggregation_stats))
    mz_stats
  }
  
  # First export: Snakemake Wildcard
  # Second export: Memory for future updates
  if (operation_mode == "Generation") {
    
    save_peptide_chunks <- peptide_chunks %>%
      mutate(Time = Sys.time()) %>%
      left_join(MS_mass_lists) %>%
      
      # Add filtering stats
      mutate(total_peptides = NA) %>%
      mutate(total_peptides = ifelse(Splice_type == "PCP", aggregation_stats$total_PCP, total_peptides)) %>%
      mutate(total_peptides = ifelse(Splice_type == "PSP", aggregation_stats$total_PSP, total_peptides)) %>%
      mutate(unique_peptides = NA) %>%
      mutate(unique_peptides = ifelse(Splice_type == "PCP", aggregation_stats$unique_PCP, unique_peptides)) %>%
      mutate(unique_peptides = ifelse(Splice_type == "PSP", aggregation_stats$unique_PSP, unique_peptides)) %>%
      left_join(mz_stats) 
    
    vroom_write(save_peptide_chunks, delim = ",", num_threads = Ncpu,
                unlist(snakemake@output[["chunk_aggregation_status"]]))
    vroom_write(save_peptide_chunks, 
                paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/", filename, ".csv"),
                delim = ",", num_threads = Ncpu, append = FALSE)
    
  } else if (operation_mode == "Update") {
    
    save_peptide_chunks <- peptide_chunks %>%
      mutate(Time = Sys.time()) %>%
      left_join(MS_mass_lists) %>%
      rbind(processed_files)  %>%
      
      # Add filtering stats
      mutate(total_peptides = NA) %>%
      mutate(total_peptides = ifelse(Splice_type == "PCP", aggregation_stats$total_PCP, total_peptides)) %>%
      mutate(total_peptides = ifelse(Splice_type == "PSP", aggregation_stats$total_PSP, total_peptides)) %>%
      mutate(unique_peptides = NA) %>%
      mutate(unique_peptides = ifelse(Splice_type == "PCP", aggregation_stats$unique_PCP, unique_peptides)) %>%
      mutate(unique_peptides = ifelse(Splice_type == "PSP", aggregation_stats$unique_PSP, unique_peptides)) %>%
      left_join(mz_stats) 
    
    vroom_write(save_peptide_chunks, delim = ",", num_threads = Ncpu,
                unlist(snakemake@output[["chunk_aggregation_status"]]))
    vroom_write(save_peptide_chunks, 
                paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/", filename, ".csv"),
                delim = ",", num_threads = Ncpu, append = FALSE)
  }
}
