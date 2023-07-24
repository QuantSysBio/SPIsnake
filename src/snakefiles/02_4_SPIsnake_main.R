### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Generate and filter peptides
#               
# input:        1. Proteome chunk.fasta
#               2. Parameters: index_length, peptide length, min_intervening_sequence length, output directory
# output:       
#               - Protein-peptide mapping .parquet
#               - Unique sequences with all filtering info (MW, RT, IC50, whether a peptide passed filters in a given mass_list) .parquet
#               - PTM peptides .parquet
#               - Peptide chunk .FASTA & .CSV for different filtering steps
#               
# author:       YH, JL, KP

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrangements))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cleaver))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # Master_table_expanded
  Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]]) %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  
  ### RT calibration
  RT_Performance_df <- fread(snakemake@input[["RT_Performance_df"]], sep = ",") %>% as_tibble()
  method = as.character(snakemake@params[["method"]])
  
  # PTMs
  max_variable_PTM = as.integer(snakemake@params[["max_variable_PTM"]])
  generate_spliced_PTMs = as.logical(snakemake@params[["generate_spliced_PTMs"]])
  PTM_chunk = as.integer(snakemake@params[["PTM_chunk"]])
  
  # Fasta input
  dir_DB_Fasta_chunks <- snakemake@params[["dir_DB_Fasta_chunks"]]
  dir_PSP_indices <- snakemake@params[["dir_PSP_indices"]]
  
  # Output dir
  dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
  dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  
  # CPUs
  active_Slurm <- ifelse(system("echo $SLURM_JOB_ID") != 0, TRUE, FALSE)
  print(active_Slurm)
  Ncpu <- snakemake@params[["max_cpus"]]
  print(paste0("Ncpu: ", Ncpu))
  # if (active_Slurm) {
  #   Ncpu <- min(parallelly::availableCores(methods = "Slurm"), Ncpu)
  # } else {
  #   Ncpu <- min(parallelly::availableCores(), Ncpu)
  # }
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  set_cpu_count(Ncpu)
  print(paste0("number of arrow CPUs: ", cpu_count()))
  retry_times = 3
  
  # DB parameters
  min_protein_length <- as.integer(snakemake@params[["min_protein_length"]])
  max_protein_length <- as.integer(snakemake@params[["max_protein_length"]])
  replace_I_with_L <- as.logical(snakemake@params[["replace_I_with_L"]])
  
  # Save into chunks according to first N letters and max netMHCpan input size
  index_length = as.integer(snakemake@params[["AA_index_length"]])
  netMHCpan_chunk = as.integer(snakemake@params[["netMHCpan_chunk"]])
  
  # netMHCpan installation path  
  netMHCpan <- as.character(snakemake@params[["netMHCpan_path"]])
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- as.logical(snakemake@params[["FASTA_outputs_unfiltered"]])
  FASTA_outputs_MW_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered"]])
  FASTA_outputs_MW_filtered_PTM <- as.logical(snakemake@params[["FASTA_outputs_MW_filtered_PTM"]])
  FASTA_outputs_MW_RT_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_filtered"]])
  FASTA_outputs_MW_RT_IC50_filtered <- as.logical(snakemake@params[["FASTA_outputs_MW_RT_IC50_filtered"]])
  
  # Whether to compress CSV output
  compress_CSV <- as.logical(snakemake@params[["compress_CSV"]])
  
  # Wildcard
  filename <- snakemake@output[[1]]
  filename_out <- filename
} else {
  ### Manual startup
  
  ### setwd("/home/yhorokh/data/Results_reports/wd/SPIsnake-3")
  # Master_table_expanded
  Master_table_expanded <- fread("results/DB_exhaustive/Master_table_expanded.csv") %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread("data/Experiment_design.csv") %>% as_tibble()
  
  ### RT calibration
  RT_Performance_df <- fread("results/RT_prediction/RT_Performance.csv", sep = ",") %>% as_tibble()
  method = as.character("achrom")
  
  # PTMs
  max_variable_PTM = 2
  generate_spliced_PTMs = FALSE
  PTM_chunk = 10000
  
  # Fasta input
  dir_DB_Fasta_chunks <- "results/DB_exhaustive/Fasta_chunks/"
  dir_PSP_indices <- "results/DB_exhaustive/PSP_indices/"
  
  # Output dir
  dir_DB_exhaustive <- "results/DB_exhaustive/"
  dir_DB_PTM_mz <- "results/DB_PTM_mz/" 
  dir_DB_out = "results/DB_out/" 
  
  # CPUs
  active_Slurm <- ifelse(system("echo $SLURM_JOB_ID") != 0, TRUE, FALSE)
  if (active_Slurm) {
    Ncpu <- min(parallelly::availableCores(methods = "Slurm") - 1)
  } else {
    Ncpu <- min(parallelly::availableCores() - 1)
  }
  print(paste0("number of CPUs: ", Ncpu))
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  retry_times = 3
  set_cpu_count(Ncpu)
  print(paste0("number of arrow CPUs: ", cpu_count()))
  
  # Save into chunks according to first N letters
  min_protein_length <- 7L
  max_protein_length <- 500L
  replace_I_with_L <- TRUE
  
  # Save into chunks according to first N letters and max netMHCpan input size
  index_length = 1L
  index_length_2 = 1
  netMHCpan_chunk = 10^6
  
  # NetMHCPan
  # netMHCpan <- "/home/yhorokh/data/Results_reports/Mouse_lymphoma_target_11.2020/bin/netMHCpan-4.1/netMHCpan"
  netMHCpan <- "/home/yhorokh/wd/netMHCpan-4.1/netMHCpan"
  
  # FASTA outputs
  FASTA_outputs_unfiltered <- TRUE
  FASTA_outputs_MW_filtered <- TRUE
  FASTA_outputs_MW_filtered_PTM <- TRUE
  FASTA_outputs_MW_RT_filtered <- TRUE
  FASTA_outputs_MW_RT_IC50_filtered <- TRUE
  
  # Whether to compress CSV output
  compress_CSV <- TRUE
  
  # Wildcard
  filename <- "results/DB_exhaustive/.5_30_trypsin_2_25_expressed_threeUTR_reversed_1.done"
  filename_out <- filename
}

# Output dirs
suppressWarnings(dir.create(paste0(dir_DB_exhaustive)))
suppressWarnings(dir.create(paste0(dir_DB_exhaustive, "/peptide_mapping")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz)))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/peptide_seqences")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/netMHCpan_input")))
suppressWarnings(dir.create(paste0(dir_DB_PTM_mz, "/netMHCpan_output")))
suppressWarnings(dir.create(paste0(dir_DB_out)))
suppressWarnings(dir.create(paste0(dir_DB_out, "/Stats")))
suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA")))
suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV")))

### Extract parameters
filename <- str_replace(filename_out, "results/DB_exhaustive/.", "results/DB_exhaustive/")
filename <- str_remove(filename, ".done") %>%
  str_split_fixed(pattern = fixed("DB_exhaustive/"), n = 2)
filename <- filename[,2]

cat(as.character(Sys.time()), " - ", "Loaded parameters", "\n")
print(filename)

params <- Master_table_expanded[Master_table_expanded$filename == filename,] %>%
  select(Proteome, enzyme_type, N_mers, Max_Interv_length, chunk, MaxE) %>%
  unique()
print(t(params))

### Input fasta + fai
{
  proteome_index <- read_parquet(paste0(dir_DB_Fasta_chunks, "/", params$Proteome, ".parquet")) %>%
    lazy_dt() %>%
    select(recno, fileno, offset, desc, seqlength, filepath, matches("maxE_", params$MaxE)) %>%
    relocate(recno, fileno, offset, desc, seqlength, filepath) %>%
    filter(if_any(matches(paste0("maxE_", params$MaxE)), ~ . == params$chunk)) %>%
    select(recno, fileno, offset, desc, seqlength, filepath) %>% 
    mutate(recno = as.integer(recno),
           seqlength = as.integer(seqlength),
           offset = as.numeric(offset)) %>%
    as.data.frame()
  cat(as.character(Sys.time()), " - ", "Reading FASTA index: Done", "\n")
  
  dat <- readAAStringSet(proteome_index)
  cat(as.character(Sys.time()), " - ", "Reading FASTA: Done", "\n")
}

# Splice type
enzyme_type <- as.character(params$enzyme_type)

# max intervening sequence length
MiSl <- as.numeric(params$Max_Interv_length)

# Nmers
if (enzyme_type == "PCP" | (grepl("cis-PSP", enzyme_type) == TRUE)) {
  # For PCP and cis-PSP corresponds to Nmer to be produced
  Nmers <- as.numeric(params$N_mers)
} else {
  # For cleaver digestions defines min and max peptide length
  Nmers <- as.numeric(unlist(str_split(params$N_mers, "_", 2)))
}

# Enzyme rules
enzymes <- c("arg-c proteinase", "asp-n endopeptidase", "bnps-skatole", "caspase1", "caspase2", "caspase3", "caspase4", "caspase5", "caspase6", "caspase7", "caspase8", "caspase9", "caspase10", "chymotrypsin-high", "chymotrypsin-low", "clostripain", "cnbr", "enterokinase", "factor xa", "formic acid", "glutamyl endopeptidase", "granzyme-b", "hydroxylamine", "iodosobenzoic acid", "lysc", "lysn", "neutrophil elastase", "ntcb", "pepsin1.3", "pepsin", "proline endopeptidase", "proteinase k", "staphylococcal peptidase i", "thermolysin", "thrombin", "trypsin")
custom_trypsin <- c("[KR](?=\\w)")

# Exclusion pattern: peptides with these letters will be omitted
exclusion_pattern <- "(U|X|O|\\*)"

### ---------------------------- (2) Inputs: MW & RT --------------------------------------
### Mass_lists
# Only user defined lists will be used 
MS_mass_lists <- list.files("data/MS_mass_lists", pattern = ".txt") %>%
  as_tibble() %>%
  dplyr::mutate(mass_list = str_remove_all(value, ".txt")) %>%
  dplyr::filter(mass_list %in% Experiment_design$Filename) %>%
  dplyr::rename(mass_list_file = value) %>%
  dplyr::mutate(AA_length = filename) 

# Load MS-1 data
MS_mass_lists_data <- list()
for (j in 1:nrow(MS_mass_lists)) {
  MS_mass_list <- MS_mass_lists$mass_list[j]
  tolerance = as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
  RT_tolerance = as.numeric(RT_Performance_df$mean_value[RT_Performance_df$dataset == MS_mass_list & RT_Performance_df$metric == "MAE"])
  MS_mass_lists_data[[j]] <- read_MW_file(file = paste0("data/MS_mass_lists/", MS_mass_lists$mass_list_file[j]), num_threads = Ncpu)
}
names(MS_mass_lists_data) <- MS_mass_lists$mass_list
rm(j)

# Masses for PTM MW filter - should contain all masses even in case of mass_list update!
mzList <- rbindlist(MS_mass_lists_data, idcol = "mzList")[,.(mzList, MW_Min, MW_Max)] %>%
  unique() %>%
  setkey(MW_Min, MW_Max, mzList)

### Retention time
# File should be named the same as the corresponding mass_list
RT_calibration_lists <- list.files("data/RT_calibration", pattern = ".csv") %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(RT_list = str_remove_all(value, ".csv")) %>%
  dplyr::filter(RT_list %in% Experiment_design$Filename) %>%
  dplyr::rename(RT_list_file = value) %>%
  dplyr::mutate(AA_length = filename) 

cat(as.character(Sys.time()), " - ", "Reading MW and RT calibration data: Done", "\n")

### ---------------------------- (3) Operation mode --------------------------------------
### Check if the peptides have been generated
chunk_params <- list(proteome = as.character(unique(params$Proteome)),
                     enzyme_type = as.character(unique(params$enzyme_type)),
                     MiSl = as.numeric(unique(params$Max_Interv_length)),
                     chunk = as.character(unique(str_remove_all(Master_table_expanded$chunk[Master_table_expanded$filename == filename],
                                                                pattern = paste0(params$Proteome, "_|.fasta")))))

if (chunk_params$enzyme_type == "cis-PSP") {
  keep_enzyme <- c("cis-PSP","PSP","PCP")
} else {
  keep_enzyme <- chunk_params$enzyme_type
}
# peptides <- list.files(paste0(dir_DB_exhaustive, "/peptide_mapping/"), pattern = ".parquet", recursive = T)
# peptides <- peptides[str_detect(peptides, str_c(paste0("enzyme=", keep_enzyme), collapse = "|"))]
# peptides <- peptides[str_detect(peptides, paste0("MiSl=", chunk_params$MiSl))]
# peptides <- peptides[str_detect(peptides, paste0("proteome=", chunk_params$proteome))]
# peptides <- peptides[str_detect(peptides, paste0("chunk=", chunk_params$chunk))]
# peptides <- peptides[str_detect(peptides, str_c(paste0("length=", seq(min(Nmers), max(Nmers))), collapse = "|"))]
# 
# keep_colnames <- c(paste0("MW:", Experiment_design$Filename),
#                    paste0("MW.exists:", Experiment_design$Filename),
#                    paste0("RT:", Experiment_design$Filename),
#                    paste0("MW.RT.exists:", Experiment_design$Filename),
#                    # paste0("Aff\\(nM\\):", Experiment_design$Filename),
#                    paste0("Predicted_binder:", Experiment_design$Filename))
# 
# if (length(peptides) == 0) {
#   cat("Generating mew peptide database\n", as.character(Sys.time()), "\n")
#   
#   ### Generate peptides
#   source("src/snakefiles/02_4_1_Generate_peptides.R")
#   
#   ### PTMs, 2D filer, IC50 prediction
#   source("src/snakefiles/02_4_2_PTM_mz_RT_matching.R")
#   
# } else {
#   DB_PTM_mz <- open_dataset(paste0(dir_DB_PTM_mz, "/peptide_seqences/")) 
#   if (!nrow(DB_PTM_mz) == 0) {
#     DB_PTM_mz <- DB_PTM_mz %>%
#       filter(enzyme %in% keep_enzyme) %>%
#       filter(MiSl == chunk_params$MiSl) %>%
#       filter(proteome == chunk_params$proteome) %>%
#       filter(chunk == as.integer(chunk_params$chunk)) %>%
#       filter(length >= min(Nmers)) %>%
#       filter(length <= max(Nmers)) 
#   }
#   
#   if (nrow(DB_PTM_mz) == 0) {
#     cat("Generating mew peptide database\n", as.character(Sys.time()), "\n")
#     
#     ### Generate peptides
#     source("src/snakefiles/02_4_1_Generate_peptides.R")
#     
#     ### PTMs, 2D filer, IC50 prediction
#     source("src/snakefiles/02_4_2_PTM_mz_RT_matching.R")
#     
#   } else {
#     cat("Found an existing peptide database\n", as.character(Sys.time()), "\n")
#     processed_DBs <- DB_PTM_mz %>%
#       head() %>%
#       collect() %>%
#       colnames() %>%
#       tidyr::expand_grid(keep_colnames) %>%
#       rename(DB_PTM_mz_col = ".") %>%
#       filter(str_starts(string = DB_PTM_mz_col, pattern = keep_colnames)) %>%
#       pull(keep_colnames) %>%
#       unique()
#     
#     keep_colnames <- keep_colnames[!keep_colnames %in% processed_DBs]
#     keep_colnames <- str_remove_all(keep_colnames, "MW:|MW.exists:|RT:|MW.RT.exists:|Aff\\(nM\\):|Predicted_binder:")
#     keep_colnames <- na.omit(unique(keep_colnames))
#     
#     ### Which mass_lists have not been processed yet
#     unprocessed_DBs <- Experiment_design %>%
#       filter(Filename %in% keep_colnames)
#     if (nrow(unprocessed_DBs) > 0) {
#       cat("Found new mass_lists to be processed:", keep_colnames, ": Starting the update!","\n", as.character(Sys.time()), "\n")
#       
#       ### Remove mass_lists that were already processed
#       Experiment_design <- unprocessed_DBs
#       MS_mass_lists <- MS_mass_lists %>% filter(mass_list %in% Experiment_design$Filename) 
#       MS_mass_lists_data <- MS_mass_lists_data[names(MS_mass_lists_data) %in% Experiment_design$Filename] 
#       RT_calibration_lists <- RT_calibration_lists %>% filter(RT_list %in% Experiment_design$Filename)  
#       
#       ### Load peptide database
#       peptides <- open_dataset(paste0(dir_DB_exhaustive, "/peptide_mapping/")) %>%
#         filter(enzyme %in% keep_enzyme) %>%
#         filter(MiSl == chunk_params$MiSl) %>%
#         filter(proteome == chunk_params$proteome) %>%
#         filter(chunk == as.integer(chunk_params$chunk))  %>%
#         filter(length >= min(Nmers)) %>%
#         filter(length <= max(Nmers)) %>%
#         select(index, peptide, length, protein, enzyme) %>%
#         collect() %>%
#         setDT()
#       setcolorder(peptides, neworder = c("index", "peptide", "protein", "enzyme", "length"))
#       cat("Loaded an existing peptide database\n", as.character(Sys.time()), "\n")
#       
#       ### PTMs, 2D filer, IC50 prediction
#       source("src/snakefiles/02_4_2_PTM_mz_RT_matching.R")
#       
#     } else {
#       cat("NO new mass_lists to be processed - terminating", "\n", as.character(Sys.time()), "\n")
#       SPIsnake_log()
#       quit(save = "no", status = 0)
#     } # no new mass_lists
#   } # peptide database exists
# }

### ---------------------------- (4) Exports ------------------------------------------------
### Generate peptides
source("src/snakefiles/02_4_1_Generate_peptides.R")

### PTMs, 2D filer, IC50 prediction
source("src/snakefiles/02_4_2_PTM_mz_RT_matching.R")

source("src/snakefiles/02_4_3_Exports.R")

### ---------------------------- (5) Log ----------------------------------------------------
SPIsnake_log()
