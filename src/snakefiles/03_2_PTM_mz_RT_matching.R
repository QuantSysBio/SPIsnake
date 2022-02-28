### ---------------------------------------------- Aggregate_PTM_mz_RT_matching ----------------------------------------------
# description:  Find a unique set of peptide sequences, compute molecular weight (MW) and do m/z and RT matching with all input mass_lists. 
#               If PTMs are required, generate them too. 
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .fst
#               - Unique peptides MW .fst
#               - Unique peptides after m/z and RT matching per experiment .fst
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
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(arrangements))

source("src/snakefiles/functions.R")
print(sessionInfo())

# create temporary directory for vroom
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

# ---------------------------- (1) Read input file and extract info ----------------------------
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
#   filename = "results/DB_PTM_mz/chunk_aggregation_status/L_8.csv"
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
#   index_length_2 = 1
#   netMHCpan_chunk = 10^5
# 
#   # RT prediction method
#   method = as.character("achrom")
# 
#   # PTMs
#   max_variable_PTM = 2
#   generate_spliced_PTMs = FALSE
#   PTM_chunk = 10000
# 
#   # bettermc::mclapply params
#   retry_times = 3
#   fst_compression = 100
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
index_length_2 = 1
netMHCpan_chunk = as.integer(snakemake@params[["netMHCpan_chunk"]])

### RT calibration data
RT_Performance_df <- vroom(snakemake@input[["RT_Performance_df"]], delim = ",", show_col_types = FALSE)

# RT prediction method
method = as.character(snakemake@params[["method"]])

# PTMs
max_variable_PTM = as.integer(snakemake@params[["max_variable_PTM"]])
generate_spliced_PTMs = as.logical(snakemake@params[["generate_spliced_PTMs"]])
PTM_chunk = as.integer(snakemake@params[["PTM_chunk"]])

### CPUs
# Ncpu = availableCores(methods = "Slurm")
Ncpu = snakemake@params[["cpus_for_R"]]
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

retry_times = 3

### fst compression level
fst_compression = as.integer(snakemake@params[["fst_compression"]])

### ---------------------------- End user variables ----------------------------
if (method == "achrom") {
  use_condaenv("R_env_reticulate")
  pyteomics <- import("pyteomics")
  
  py_run_string("
import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None
")}

# Proteinogenic AAs
AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

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

### ---------------------------- (2) Operation mode --------------------------------------
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_for_NetMHCpan/")))
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

### ---------------------------- (3) Define file processing --------------------------------------
# Select all files with the same {AA-index}{length}
peptide_chunks <- list.files(paste0(dir_DB_exhaustive, "/peptide_seqences"), pattern = ".fst", full.names = T) %>%
  as_tibble() %>%
  mutate(file = str_remove_all(value, ".fst")) %>%
  mutate(file = str_split_fixed(file, "/peptide_seqences/", n = 2)[,2]) %>%
  mutate(Splice_type = str_split_fixed(file, "_", n = 2)[,1]) %>%
  mutate(file = str_sub(file, start = nchar(Splice_type) + 2)) %>%
  mutate(AA_length = str_c(str_split_fixed(file, "_", n = 3)[,1],
                           str_split_fixed(file, "_", n = 3)[,2], sep = "_")) %>%
  filter(AA_length == filename) %>%
  # add proteome
  mutate(Proteome = str_extract(file, str_c(unique(Master_table_expanded$Proteome), collapse = "|"))) %>%
  mutate(Proteome = str_remove_all(Proteome, regex("_[\\d+]+_[\\d+]+.fasta$"))) %>%
  # add PTM
  mutate(AA = paste0(str_split_fixed(AA_length, "_", 2)[,1],"_")) %>%
  mutate(filename = str_split_fixed(file, pattern = AA, 2)[,2]) %>%
  left_join(select(Master_table_expanded, filename, PTMs))

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

### FOR TESTING !
# if ("SwissProt_UP000005640" %in% peptide_chunks$Proteome) {
#   peptide_chunks <- peptide_chunks[peptide_chunks$Splice_type == "PCP",]
#   peptide_chunks <- peptide_chunks[peptide_chunks$Proteome == "proteome_expressed_gencode",]
# Proteome_i = "proteome_expressed_gencode"
# enzyme_type = "PCP"
# PTM = "5-defaults"
# peptide_chunks <- peptide_chunks %>%
#   mutate(PTMs = ifelse(PTMs == "5-defaults", NA, PTMs))
# }

# End the script if no files to be updated
if (nrow(peptide_chunks) == 0) {
  print("Looks like there are no new peptides to be processed")
  save_peptide_chunks <- peptide_chunks %>%
    rbind(processed_files) %>%
    mutate(Time = Sys.time()) 
  
  vroom_write(save_peptide_chunks, delim = ",", num_threads = Ncpu, append = TRUE,
              unlist(snakemake@output[["chunk_aggregation_status"]]))
} else {
  Proteomes <- peptide_chunks$Proteome %>% unique()
  enzyme_types <- peptide_chunks$Splice_type %>% unique()
  mz_stats <- tibble()
  
  # Load MS-1 data
  MS_mass_lists_data <- list()
  for (j in 1:nrow(MS_mass_lists)) {
    MS_mass_list <- MS_mass_lists$mass_list[j]
    tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
    RT_tolerance = as.numeric(RT_Performance_df$mean_value[RT_Performance_df$dataset == MS_mass_list & RT_Performance_df$metric == "MAE"])
    MS_mass_lists_data[[j]] <- read_MW_file(file = paste0("data/MS_mass_lists/", MS_mass_lists$mass_list_file[j]), num_threads = Ncpu)
  }
  names(MS_mass_lists_data) <- MS_mass_lists$mass_list
  
  ### ---------------------------- Start filtering  --------------------------------------
  filename_exports <- c()
  for (Proteome_i in Proteomes) {
    peptide_chunks_tmp <- peptide_chunks[peptide_chunks$Proteome == Proteome_i,]
    enzyme_types <- peptide_chunks_tmp$Splice_type %>% unique()
    
    for (enzyme_type in enzyme_types) {
      peptides <- peptide_chunks_tmp$value[peptide_chunks_tmp$Splice_type == enzyme_type]
      names(peptides) <- peptide_chunks_tmp$file[peptide_chunks_tmp$Splice_type == enzyme_type]
      
      ### Load peptide sequences
      input <- peptides %>%
        bettermc::mclapply(mc.cores = Ncpu, mc.retry = retry_times, mc.cleanup=T, mc.preschedule=T, 
                           FUN = read_fst, as.data.table = TRUE, columns = "peptide") %>%
        rbindlist() %>%
        setkey(cols = "peptide") %>%
        unique()
      cat("Reading in", Proteome_i, enzyme_type, "sequences: Done\n",
          as.character(Sys.time()), "\n")

      ### ---------------------------- (4) Compute MW --------------------------------------
      ### Define chunking strategy
      chunk_size = PTM_chunk
      n_chunks = ceiling((nrow(input) / chunk_size))
      print(paste0("Creating ", n_chunks, " chunks"))       

      ### Compute MW block-wise
      input[, chunks := rep(1:n_chunks, each=chunk_size, length.out = nrow(input))] %>%
        .[, MW := computeMZ_biostrings(peptide), by = chunks]
      print(paste0("Created ", n_chunks, " chunks"))  
      
      cat("MW", enzyme_type, ": Done\n",
          as.character(Sys.time()), "\n")
      ### ---------------------------- (5) PTM generation --------------------------------------
      # PTMs to be generated
      PTM_list <- peptide_chunks_tmp$PTMs[peptide_chunks_tmp$Splice_type == enzyme_type] %>%
        str_split(pattern = "\\|") %>%
        unlist() %>%
        str_squish()
      PTM_list <- PTM_list[!(PTM_list == "")] %>% unique() %>% na.omit() %>% as.character()
      PTM_list <- PTM_list[!PTM_list == "OpenSearch"]
      
      if (generate_spliced_PTMs == FALSE) {
        ignore_generate_spliced_PTMs = !(enzyme_type == "PSP" | enzyme_type == "cis-PSP")
      } else {
        ignore_generate_spliced_PTMs = TRUE
      }
      
      if (length(PTM_list) > 0 & ignore_generate_spliced_PTMs) {
        mzList <- rbindlist(MS_mass_lists_data, idcol = "mzList")[,.(mzList, MW_Min, MW_Max)] %>%
          unique() %>%
          setkey(MW_Min, MW_Max, mzList)
        
        for (PTM in PTM_list) {
          # Load PTM table
          mods <- vroom(paste0("data/modifications/", PTM, ".csv"), show_col_types = F, num_threads = Ncpu) 
          
          # Prepare peptides
          counter_fst = 1
          for (input_i in 1:n_chunks) {
            print(paste("Generating PTMs for:", PTM, input_i))
            
            # Additional split to reduce RAM usage
            PTMcombinations <- input[chunks == input_i, .(peptide, MW)]
            PTMcombinations <- PTMcombinations[, PTMcombinations := parallel::mcmapply(getPTMcombinations_fast_vec, 
                                                                                       peptide, MW, max_variable_PTM, 
                                                                                       list(mods), 
                                                                                       SIMPLIFY = T, 
                                                                                       mc.cores = Ncpu, 
                                                                                       mc.preschedule = T, 
                                                                                       USE.NAMES = F, 
                                                                                       mc.cleanup = T)]
            print("got PTM combinations")
            
            PTMcombinations <- rbindlist(PTMcombinations$PTMcombinations, use.names = F) %>%
              .[!is.na(ids)] %>%
              .[, id := rowid(peptide, ids)] %>%
              .[, c("MW_Min", "MW_Max") := MW]
            print("rbind combinations")
            
            # If there have been PTMs generated
            if (!(nrow(PTMcombinations) == 1 & anyNA(PTMcombinations$ids))) {
              PTM_stats_all <- PTMcombinations[, .(All_PTM = .N,
                                                   Unique_PTM = uniqueN(ids)), keyby = .(peptide)]
              
              PTM_MW_out <- foverlaps(x = PTMcombinations, y = mzList, 
                                      by.x=c("MW_Min","MW_Max"),
                                      by.y=c("MW_Min", "MW_Max"), nomatch=0L) %>%
                .[, c("i.MW_Min", "i.MW_Max", "MW_Min", "MW_Max"):=NULL] %>%
                unique() %>%
                .[, id := NULL]
              print("filtered by MW")
              
              # Save modified peptides
              counter_fst = counter_fst + 1
              PTM_MW_out %>%
                vroom_write(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/", enzyme_type, "_", filename, "_", Proteome_i, "_", PTM, ".csv.gz"), 
                            delim = ",", append = T, num_threads = Ncpu)
              
              # Save stats
              PTM_MW_out %>%
                .[, .(All_MW_filtered_PTM = .N,
                      Unique_MW_filtered_PTM = uniqueN(ids)), keyby = .(peptide, mzList)] %>%
                .[PTM_stats_all, on = .(peptide), allow.cartesian = TRUE] %>%
                relocate(peptide, All_PTM, Unique_PTM, All_MW_filtered_PTM, Unique_MW_filtered_PTM, mzList) %>%
                vroom_write(paste0(
                  dir_DB_PTM_mz, "/stats_PTM/", enzyme_type, "_", filename, "_", Proteome_i, "_", PTM, ".csv.gz"), 
                  delim = ",", append = T, num_threads = Ncpu)
            } 
            rm(PTMcombinations, PTM_MW_out, PTM_stats_all)
          } # end input_i
        } # end PTM
      } # end PTM generation
      ### ---------------------------- (6) MW & RT matching --------------------------------------
      for (j in 1:nrow(MS_mass_lists)) {
        print(MS_mass_lists$mass_list[j])
        MS_mass_list <- MS_mass_lists$mass_list[j]
        
        # Filter tolerances
        tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
        RT_tolerance = as.numeric(RT_Performance_df$mean_value[RT_Performance_df$dataset == MS_mass_list & RT_Performance_df$metric == "MAE"])
        
        # Calibration input
        mzList = MS_mass_lists_data[[j]]
        
        # MW filter block-wise
        tmp <- input[MW %inrange% mzList[,c("MW_Min", "MW_Max")], .(peptide, MW), by = chunks]
        
        # Save stats
        mz_stats_ij = tibble(Proteome = Proteome_i,
                             Splice_type = enzyme_type,
                             mass_list = MS_mass_list,
                             unique_peptides = nrow(input),
                             mz_matched_peptides = nrow(tmp))
        cat("MW filter: Done\n",
            as.character(Sys.time()), "\n")

        ### ----------------------------- Predict RT for m/z matched peptides -----------------------------
        if (method == "AutoRT") {
          # Save AutoRT input
          tmp <- tmp %>%
            bettermc::mclapply(mc.cores = 1, mc.retry = retry_times, mc.cleanup=T, mc.preschedule=T, FUN = function(x){
                  # Save peptides for AutoRT prediction
                  x %>%
                    as_tibble() %>%
                    rename(x=peptide) %>%
                    select(x) %>%
                    vroom_write(file = paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", enzyme_type, "_", filename, "_", Proteome_i, "_", MS_mass_list, ".tsv"), 
                                num_threads = Ncpu, append = F, delim = "\t")
                  
                  # AutoRT predict with pre-trained model
                  system(command = paste("python bin/AutoRT/autort.py predict --test", 
                                         paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", enzyme_type, "_", filename, "_", Proteome_i, "_", MS_mass_list, ".tsv"),
                                         "-s", paste0("results/RT_prediction/AutoRT_models/", MS_mass_list, "/model.json"),
                                         "-o", paste0("results/RT_prediction/peptide_RT/", enzyme_type, "_", filename, "_", MS_mass_list)), 
                         intern = T)
                  
                  # Import AutoRT
                  RT_pred <- vroom(paste0("results/RT_prediction/peptide_RT/", enzyme_type, "_", filename, "_", Proteome_i, "_", MS_mass_list, "/test.tsv"), show_col_types = FALSE)
                  x$RT_pred <- RT_pred$y_pred
                  return(x)
                })
        } else if (method == "achrom") {
          # Save unique peptides
          tmp %>%
            as_tibble() %>%
            select(-chunks) %>%
            vroom_write(file = paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", enzyme_type, "_", filename, "_", Proteome_i, "_", MS_mass_list, ".tsv"), 
                        num_threads = Ncpu, append = T, delim = "\t")
          
          # Load parameters
          RCs <- readRDS(paste0("results/RT_prediction/RT_models/", MS_mass_list, "/achrom_RCs.rds"))
          
          # Use Leu values for Ile in case no Ile was provided in calibration
          if (is.null(RCs$aa$I)) {
            RCs$aa$I <- RCs$aa$L
          }
          
          # Exclusion pattern: peptides containing these letters will be omitted
          exclusion_pattern <- AA[!AA %in% names(RCs$aa)] %>% str_c(collapse = "|")
          if (!exclusion_pattern == "") {
            tmp <- tmp[!str_detect(peptide, exclusion_pattern), .(peptide, MW), by = chunks]
          }
          
          ### Predict RT
        not_empty_MS_mass_list <- nrow(tmp) > 0
          
          if (not_empty_MS_mass_list) {
            py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
")          
            tmp <- tmp[, 
                       RT_pred := split(peptide, factor(ceiling(seq_along(peptide)/ceiling(length(peptide)/Ncpu)))) %>%
                         bettermc::mclapply(mc.cores = Ncpu, mc.retry = retry_times, mc.cleanup=T, mc.preschedule=T, 
                                            FUN = function(x){
                                              if (length(x) == 1) {
                                                unique(py_calls$achrom_calculate_RT(rep(x, 2), RCs, r_to_py(FALSE)))
                                              } else {
                                                py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE))
                                              }
                                            }) %>%
                         unlist()]
          }
        }
        cat("RT prediction: Done\n",
            as.character(Sys.time()), "\n")
        
        ### 2D filter: MW & RT
        not_empty_MS_mass_list <- nrow(tmp) > 0
        if (not_empty_MS_mass_list) {

          tmp <- mzList[tmp, on=.(MW_Min <= MW, MW_Max >= MW,
                                   RT_Min <= RT_pred, RT_Max >= RT_pred), nomatch=0,
                         .(peptide, MW, RT_pred)]
        }
        cat("2D MW/RT filter: Done\n",
            as.character(Sys.time()), "\n")
        
        # Update stats after RT filter
        mz_stats_ij$mz_RT_matched_peptides = nrow(tmp)
        mz_stats = rbind(mz_stats, mz_stats_ij)
        
        ### Save filtered peptides
        if (nrow(tmp) > 0) {
          filename_exports <- c(filename_exports,
                                paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", enzyme_type, "_", 
                                       filename, "_", Proteome_i, "_", MS_mass_lists$mass_list[j],".fst"))
          threads_fst(nr_of_threads = Ncpu)
          write_fst(tmp, path = last(filename_exports), compress = fst_compression)
        }
        rm(tmp)
      } # End 2D filter
    } # End enzyme type
  } # End proteome
  rm(input)
  
  ### ----------------------------- Prepare NetMHCPan inputs ----------------------------- 
  # filename_exports <- list.files("results/DB_PTM_mz/unique_peptides_mz_RT_matched/", pattern = "L_8", full.names = T)
  filename_exports <- filename_exports[file.exists(filename_exports)]
  names(filename_exports) <- str_split_fixed(filename_exports, "/unique_peptides_mz_RT_matched/", 2)[,2] %>%
    str_remove(".fst")
  
  # Define IC50 aggregation
  IC50_aggregation_table <- Experiment_design %>%
    filter(!(is.na(`MHC-I_alleles`) | is.na(Affinity_threshold))) %>%
    tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]", ) %>%
    mutate(`MHC-I_alleles` = str_squish(`MHC-I_alleles`),
           Affinity_threshold = str_squish(Affinity_threshold)) 
  IC50_aggregation_table
  
  alleles <- IC50_aggregation_table$`MHC-I_alleles` %>% unique() %>% na.omit()
  for (allele in alleles) {
    IC50_aggregation <- IC50_aggregation_table$Filename[IC50_aggregation_table$`MHC-I_alleles` == allele] %>%
      str_c(collapse = "|")
    keep <- str_ends(names(filename_exports), IC50_aggregation)
    
    if (TRUE %in% keep) {
      # Export all unique peptides for a given allele
      peptide_exports <- filename_exports[keep] %>%
        bettermc::mclapply(mc.cores = Ncpu, mc.retry = retry_times, mc.cleanup=T, mc.preschedule=T, 
                           FUN = read_fst, as.data.table = TRUE, columns = "peptide") %>%
        rbindlist() %>%
        setkey("peptide")
      
        peptide_exports[, netMHCpan_split := ceiling(seq_along(peptide)/netMHCpan_chunk)] %>%
          split(by = "netMHCpan_split", drop = T, keep.by = T) %>%
          bettermc::mclapply(mc.cores = 2, mc.cleanup=T, mc.preschedule=F, FUN = function(x){
              netMHCpan_split <- x$netMHCpan_split[[1]]
              
              unique(x)[, .(peptide)] %>% 
                vroom_write(file = paste0(dir_DB_PTM_mz, "/unique_peptides_for_NetMHCpan/", 
                                          filename, "_", allele, "_ch_", netMHCpan_split ,".tsv"),
                            delim = "\t", num_threads = Ncpu, append = FALSE, col_names = FALSE)            
              
          }) %>%
          capture.output()
    }
  } 
  ### --------------------------------------(8) Chunk aggregation status (Snakemake output) --------------------------------------
  # First export: Snakemake Wildcard
  # Second export: Memory for future updates
  if (operation_mode == "Generation") {
    
    save_peptide_chunks <- peptide_chunks %>%
      mutate(Time = Sys.time()) %>%
      left_join(MS_mass_lists) %>%
      left_join(mz_stats) 
    
  } else if (operation_mode == "Update") {
    
    save_peptide_chunks <- peptide_chunks %>%
      mutate(Time = Sys.time()) %>%
      left_join(MS_mass_lists) %>%
      left_join(mz_stats) %>%
      rbind(processed_files) 
  }
  vroom_write(save_peptide_chunks, delim = ",", num_threads = Ncpu,
              unlist(snakemake@output[["chunk_aggregation_status"]]))
  vroom_write(save_peptide_chunks, 
              paste0(dir_DB_PTM_mz, "/chunk_aggregation_memory/", filename, ".csv"),
              delim = ",", num_threads = Ncpu, append = FALSE)
  
} # End peptide processing

### --------------------------------------(9) Cluster termination & stats --------------------------------------
print("----- removing cluster -----")
print(cl)
parallel::stopCluster(cl)

# rm(rf=list())
gc()

print(Sys.time())
print(paste("Finished MW & RT filtering"))

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

print("----- garbage collection -----")
gc() %>%
  print()
