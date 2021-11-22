### ---------------------------------------------- Aggregate_PTM_mz_RT_matching ----------------------------------------------
# description:  Find a unique set of peptide sequences, compute molecular weight (MW) and do m/z and RT matching with all input mass_lists. 
#               If PTMs are required, generate them too. 
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .csv.gz
#               - Unique peptides MW .csv.gz
#               - Unique peptides after m/z and RT matching per experiment .csv.gz
#               - Peptide-mass_list matching
#               - Peptide filtering stats .csv
#               
# author:       YH, JL, KP

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(arrangements))

source("src/snakefiles/functions.R")

### ---------------------------- (1) Read input file and extract info ----------------------------
# {
#   ### Manual setup
#   # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
# 
#   Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv")
#   Peptide_aggregation_table <- vroom("results/DB_PTM_mz/Peptide_aggregation_table.csv", delim = ",")
#   Experiment_design <- vroom("data/Experiment_design.csv", delim = ",")
#   dir_DB_exhaustive = "results/DB_exhaustive"
#   dir_DB_PTM_mz = "results/DB_PTM_mz"
#   suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory")))
# 
#   # Wildcard
#   filename = "results/DB_PTM_mz/chunk_aggregation_status/V_11.csv"
#   filename <- filename %>%
#     str_split_fixed(pattern = fixed("chunk_aggregation_status/"), n = 2)
#   filename <- filename[,2] %>%
#     str_remove(pattern = ".csv")
#   print(filename)
# 
#   # Calibration
#   MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
#     as_tibble() %>%
#     mutate(file = str_remove_all(value, ".txt"))
#   RT_Performance_df <- vroom("results/RT_prediction/RT_Performance.csv", delim = ",", show_col_types = FALSE)
# 
#   ### CPUs
#   Ncpu = availableCores()
#   cl <- parallel::makeForkCluster(Ncpu)
#   cl <- parallelly::autoStopCluster(cl)
#   setDTthreads(Ncpu)
# 
#   # Save into chunks according to first N letters
#   index_length = 1
#   netMHCpan_chunk = 10^5
# 
#   # RT prediction method
#   method = as.character("achrom")
# 
#   # PTMs
#   max_variable_PTM = 2
#   generate_spliced_PTMs = FALSE
# }

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

# Save into chunks according to first N letters and max netMHCpan input size
index_length = as.integer(snakemake@params[["AA_index_length"]])
netMHCpan_chunk = as.integer(snakemake@params[["netMHCpan_chunk"]])

### RT calibration data
RT_Performance_df <- vroom(snakemake@input[["RT_Performance_df"]], delim = ",", show_col_types = FALSE)

# RT prediction method
method = as.character(snakemake@params[["method"]])

# PTMs
max_variable_PTM = as.integer(snakemake@params[["max_variable_PTM"]])
generate_spliced_PTMs = as.logical(snakemake@params[["generate_spliced_PTMs"]])

### ---------------------------- End user variables ----------------------------
# Proteinogenic AAs
AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

### Mass_lists
# Only user defined lists will be used 
MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
  as_tibble() %>%
  mutate(mass_list = str_remove_all(value, ".txt")) %>%
  filter(mass_list %in% Experiment_design$Filename) %>%
  rename(mass_list_file = value) %>%
  mutate(AA_length = filename) 

# File should be named the same as the corresponding mass_list
RT_calibration_lists <- list.files("data/RT_calibration", pattern = ".csv") %>%
  as_tibble() %>%
  mutate(RT_list = str_remove_all(value, ".csv")) %>%
  filter(RT_list %in% Experiment_design$Filename) %>%
  rename(RT_list_file = value) %>%
  mutate(AA_length = filename) 

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (2) Operation mode --------------------------------------
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/stats_PTM/")))

# Check if there exist previous outputs to be updated:
processed_files <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory"), pattern = paste0(filename, ".csv"))

operation_mode = ifelse(length(processed_files) > 0, "Update", "Generation")
print(Sys.time())
print(paste0("Mode: ", operation_mode))

if (operation_mode == "Update") {
  print("Updating file list to be processed")
  processed_files <- vroom(file = paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/", processed_files), delim = ",",  show_col_types = FALSE)
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
  # peptides <- peptides[1:60]
  
  # PCP
  if (length(peptides[str_detect(peptides, "/PCP_")]) > 0) {
    PCP_list <- peptides[str_detect(peptides, "/PCP_")]  %>%
      mclapply(FUN = vroom, delim = ",", mc.cores = Ncpu, show_col_types = FALSE) 
    print(Sys.time())
    print("Reading in PCP sequences: Done")
  }
  # PSP
  if (length(peptides[str_detect(peptides, "/PSP_")]) > 0) {
    PSP_list <- peptides[str_detect(peptides, "/PSP_")]  %>%
      mclapply(FUN = vroom, delim = ",", mc.cores = Ncpu, show_col_types = FALSE)
    print(Sys.time())
    print("Reading in PSP sequences: Done")
  }
  
  ### ---------------------------- (4) Compute MW --------------------------------------
  ### And check for peptide uniqueness
  pep_types <- c()
  if (length(peptides[str_detect(peptides, "/PCP_")]) > 0) {
    PCP <- PCP_list %>%
      rbindlist()
    PCP$index <- str_sub(PCP$peptide, index_length + 1, index_length + 1)
    PCP <- setorder(PCP, index)
    
    PCP <- PCP %>%
      split(by = "index", drop = T) %>%
      mclapply(mc.cores = Ncpu, FUN = function(x){
        x %>%
          lazy_dt() %>%
          unique() %>%
          mutate(MW=computeMZ_biostrings(peptide)) %>%
          as.data.table()
      }) 
    pep_types <- c(pep_types, "PCP")
    print(Sys.time())
    print("MW PCP: Done")
  }
  
  if (length(peptides[str_detect(peptides, "/PSP_")]) > 0) {
    PSP <- PSP_list %>%
      rbindlist() 
    PSP$index <- str_sub(PSP$peptide, index_length + 1, index_length + 1)
    PSP <- setorder(PSP, index)
    
    PSP <- PSP %>%
      split(by = "index", drop = T) %>%
      mclapply(mc.cores = Ncpu, FUN = function(x){
        x %>%
          lazy_dt() %>%
          unique() %>%
          mutate(MW=computeMZ_biostrings(peptide)) %>%
          as.data.table()
      }) 
    pep_types <- c(pep_types, "PSP")
    print(Sys.time())
    print("MW PSP: Done")
  }
  
  ### ---------------------------- (5) MW & RT matching --------------------------------------
  if (method == "achrom") {
    use_condaenv("R_env_reticulate")
    pyteomics <- import("pyteomics")
    
    py_run_string("
import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None
")  
  }
  
  mz_nomod <- list()
  mz_stats <- tibble()
  for (i in pep_types) {
    input <- get(i)
    
    for (j in 1:nrow(MS_mass_lists)) {
      print(MS_mass_lists$mass_list[j])
      MS_mass_list <- MS_mass_lists$mass_list[j]
      
      # Filter tolerances
      tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
      RT_tolerance = as.numeric(RT_Performance_df$mean_value[RT_Performance_df$dataset == MS_mass_list & RT_Performance_df$metric == "MAE"])
      
      # Calibration input
      mzList = read_MW_file(file = paste0("data/MS_mass_lists/", MS_mass_lists$mass_list_file[j]), num_threads = Ncpu)
      
      # Which peptide sequences pass the MW filter
      mz_nomod[[i]][[MS_mass_list]] <- input %>%
        mclapply(mc.cores = Ncpu, FUN = function(x){
          x[MW %inrange% mzList[,c("MW_Min", "MW_Max")]]
        }) %>%
        rbindlist()
      print(Sys.time())
      print("MW filter: done")
      
      # Save stats
      mz_stats_ij = tibble(Splice_type=i,
                           mass_list = MS_mass_list,
                           mz_matched_peptides = nrow(mz_nomod[[i]][[MS_mass_list]]))
      
      ### ----------------------------- Predict RT for m/z matched peptides -----------------------------
      if (method == "AutoRT") {
        # Save AutoRT input
        mz_nomod[[i]][[MS_mass_list]]  %>%
          as_tibble() %>%
          rename(x=peptide) %>%
          select(x) %>%
          vroom_write(file = paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", i, "_", filename, "_", MS_mass_list, ".tsv"), 
                      num_threads = Ncpu, append = F, delim = "\t")
        
        # AutoRT predict with pre-trained model
        system(command = paste("python bin/AutoRT/autort.py predict --test", 
                               paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", i, "_", filename, "_", MS_mass_list, ".tsv"),
                               "-s", paste0("results/RT_prediction/AutoRT_models/", MS_mass_list, "/model.json"),
                               "-o", paste0("results/RT_prediction/peptide_RT/", i, "_", filename, "_", MS_mass_list)), 
               intern = T)
        
        # Import AutoRT
        RT_pred <- vroom(paste0("results/RT_prediction/peptide_RT/", i, "_", filename, "_", MS_mass_list, "/test.tsv"), show_col_types = FALSE)
        mz_nomod[[i]][[MS_mass_list]]$RT_pred <- RT_pred$y_pred
        
      } else if (method == "achrom") {
        mz_nomod[[i]][[MS_mass_list]]  %>%
          as_tibble() %>%
          vroom_write(file = paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", i, "_", filename, "_", MS_mass_list, ".tsv"), 
                      num_threads = Ncpu, append = F, delim = "\t")
        
        # Load parameters
        RCs <- readRDS(paste0("results/RT_prediction/RT_models/", MS_mass_list, "/achrom_RCs.rds"))
        
        # Use Leu values for Ile in case no Ile was provided in calibration
        if (is.null(RCs$aa$I)) {
          RCs$aa$I <- RCs$aa$L
        }
        
        # Exclusion pattern: peptides containing these letters will be omitted
        exclusion_pattern <- AA[!AA %in% names(RCs$aa)] %>% str_c(collapse = "|")
        if (!exclusion_pattern == "") {
          mz_nomod[[i]][[MS_mass_list]] <- mz_nomod[[i]][[MS_mass_list]] %>%
            lazy_dt() %>%
            filter(!str_detect(peptide, exclusion_pattern)) %>%
            as.data.table()
        }
        
        if (!is.na(nrow(mz_nomod[[i]][[MS_mass_list]]))) {
          ### Predict
          py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
")
          mz_nomod[[i]][[MS_mass_list]]$RT_pred <- mz_nomod[[i]][[MS_mass_list]] %>%
            split(by = c("index"), drop = T) %>%
            mclapply(mc.cores = Ncpu, FUN = function(x){
              pep = x %>%
                lazy_dt() %>%
                select(-index) %>%
                pull(peptide) 
              if (length(pep) == 1) {
                pep = c(pep, pep)
                x = data.table(RT = as.numeric(py_calls$achrom_calculate_RT(pep, RCs, r_to_py(FALSE)))) %>%
                  unique()
              } else {
                x = data.table(RT = as.numeric(py_calls$achrom_calculate_RT(pep, RCs, r_to_py(FALSE))))
              }
              return(x)
            }) %>%
            rbindlist() %>%
            pull(RT)
        }
      }
      print(Sys.time())
      print("RT prediction: Done")
      
      ### 2D filter: MW & RT
      if (!is.na(nrow(mz_nomod[[i]][[MS_mass_list]]))) {
        
        mz_nomod[[i]][[MS_mass_list]] <- mz_nomod[[i]][[MS_mass_list]] %>%
          split(by = c("index"), drop = T) %>%
          mclapply(mc.cores = Ncpu, FUN = function(x){
            x = x %>%
              lazy_dt() %>%
              select(-index) %>%
              as.data.table()
            x = x[MW %inrange% mzList[,c("MW_Min", "MW_Max")] & RT_pred %inrange% mzList[,c("RT_Min", "RT_Max")]]
            return(x)
          }) %>%
          rbindlist()
      }
      print(Sys.time())
      print("2D MW/RT filter: Done")
      
      # Update stats after RT filter
      mz_stats_ij$mz_RT_matched_peptides = nrow(mz_nomod[[i]][[MS_mass_list]])
      mz_stats = rbind(mz_stats, mz_stats_ij)
    }
  }
  
  ### ---------------------------- (6) Generate PTMs --------------------------------------
  # Define PTM generation jobs
  {
    tmp <- Master_table_expanded %>%
      mutate(Proteome_PTM = paste(PTMs, Min_Interv_length, chunk, sep = "_")) %>%
      select(-Splice_type) 
    
    Master_table_expanded_PTM <- peptide_chunks %>%
      mutate(Proteome_PTM = str_split_fixed(file, pattern = AA_length, 2)[,2]) %>%
      mutate(Proteome_PTM = ifelse(str_starts(Proteome_PTM, "_"), str_sub(Proteome_PTM, 2), Proteome_PTM)) %>%
      left_join(tmp)
    Master_table_expanded_PTM[1,] %>% t()
    
    PTM_list <- Master_table_expanded_PTM$PTMs %>%
      unique() %>%
      na.omit()
    print(Sys.time())
    print(paste("Starting PTM generation for the following PTMs:", str_flatten(PTM_list, collapse = ", ")))
  }
  
  if (length(PTM_list) > 0) {
    for (PTM in PTM_list) {
      # Load PTM table
      df <- Master_table_expanded_PTM %>%
        filter(PTMs == PTM)
      mods <- vroom(paste0("data/modifications/", PTM, ".csv"), show_col_types = F) 
      
      # Which peptide types undergo PTM generation
      pep_types_PTM <- c()
      if (length(PCP_list[unique(df$file)]) > 0) {
        pep_types_PTM <- c(pep_types_PTM, "PCP")
      }
      
      # Whether to generate PTMs for spliced peptides
      if (generate_spliced_PTMs == TRUE) {
        if (length(PSP_list[unique(df$file)]) > 0) {
          pep_types_PTM <- c(pep_types_PTM, "PSP")
        }
      }
      
      for (i in pep_types_PTM) {
        
        # Prepare all relevant proteomes
        input <- get(paste0(i, "_list"))[unique(df$file)] %>%
          rbindlist()
        input$index <- str_sub(input$peptide, index_length + 1, index_length + 1)
        input <- input %>%
          split(by = "index", drop = T) %>%
          mclapply(mc.cores = Ncpu, FUN = function(x){
            x %>%
              lazy_dt() %>%
              unique() %>%
              mutate(MW=computeMZ_biostrings(peptide)) %>%
              as.data.table()
          }) %>%
          rbindlist() %>%
          split(by = "peptide", drop = T)
        
        for (j in 1:nrow(MS_mass_lists)) {
          print(MS_mass_lists$mass_list[j])
          MS_mass_list <- MS_mass_lists$mass_list[j]
          tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
          
          # MW input
          mzList = read_MW_file(file = paste0("data/MS_mass_lists/", MS_mass_lists$mass_list_file[j]), num_threads = Ncpu)
          
          # PTM generation and MW filtering per peptide
          tmp = input %>%
            mclapply(mc.cores = Ncpu, FUN = function(x){
              PTMcombinations = getPTMcombinations_fast(c(x$peptide, x$MW), NmaxMod = max_variable_PTM, mods_input = mods)
              y = PTMcombinations[MW %inrange% mzList[,c("MW_Min", "MW_Max")]]
              PTM_pep_stats = tibble(peptide = x$peptide,
                                     All_PTM = nrow(PTMcombinations),
                                     MW_filtered_PTM = nrow(y))
              return(list(pep = y,
                          PTM_pep_stats = PTM_pep_stats))
            })
          
          # Save modified peptides
          rbindlist(lapply(tmp, `[[`, 1)) %>%
            vroom_write(file = paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/", i, "_", filename, "_", PTM, "_", MS_mass_list, ".csv.gz"),
                        delim = ",", num_threads = Ncpu, append = FALSE)
          
          # Save PTM generation stats
          rbindlist(lapply(tmp, `[[`, 2)) %>%
            vroom_write(file = paste0(dir_DB_PTM_mz, "/stats_PTM/", i, "_", filename, "_", PTM, "_", MS_mass_list, ".csv.gz"),
                        delim = ",", num_threads = Ncpu, append = FALSE)
          rm(tmp)
        }
      }
    }
  }
  
  ### ---------------------------- (7) Peptide sequence export --------------------------------------
  # Unmodified m/z & RT matched sequences
  # No compression is done to facilitate read-in by NetMHCPan
    print(Sys.time())
    print(paste("Starting peptide sequence export"))
  
  for (i in seq_along(mz_nomod)) {
    for (j in seq_along(mz_nomod[[i]])) {
      mz_nomod[[i]][[j]] %>%
        lazy_dt() %>%
        mutate(netMHCpan_split = ceiling(seq_along(peptide)/netMHCpan_chunk)) %>%
        group_by(netMHCpan_split) %>%
        group_walk(~ vroom_write(.x, 
                                 file = paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", names(mz_nomod)[i], "_", filename, "_", names(mz_nomod[[i]])[j], "_ch_", .y$netMHCpan_split ,".tsv"),
                                 delim = "\t", num_threads = Ncpu, append = TRUE))
    }
  }
  
  ### --------------------------------------(8) Chunk aggregation status (Snakemake output) --------------------------------------
  ### Record aggregation stats
    print(Sys.time())
    print(paste("Starting aggregation stats export"))
  
  if (length(peptides[str_detect(peptides, "/PSP_")]) > 0) {
    aggregation_stats <- tibble(total_PCP = sum(unlist(lapply(PCP_list, nrow))), 
                                unique_PCP = sum(unlist(lapply(PCP, nrow))),
                                total_PSP = sum(unlist(lapply(PSP_list, nrow))),
                                unique_PSP = sum(unlist(lapply(PSP, nrow))))
  } else {
    aggregation_stats <- tibble(total_PCP = sum(unlist(lapply(PCP_list, nrow))), 
                                unique_PCP = sum(unlist(lapply(PCP, nrow))),
                                total_PSP = 0,
                                unique_PSP = 0)
  }
  print(t(aggregation_stats))
  
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

print(cl)
rm(list = "cl")
gc()

print(Sys.time())
print(paste("Finished MW & RT filtering"))

