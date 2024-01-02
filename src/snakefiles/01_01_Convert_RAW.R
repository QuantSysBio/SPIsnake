### ---------------------------------------------- Convert RAW to MGF  ----------------------------------------------
# description:  use ThermoRawFileParser to create MGF files
#               
# input:        1. Experiment design table specifies which files to process
#               2. ThermoRawFileParser.exe
#               3. PEAKS search outputs
#
# output:       
#               For every specified file:
#               - an .mgf file
#               - a text file with observed masses and retention times
#               - a text file with detected peptide sequences and retention times
#               
# author:       YH, JAC

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

pyteomics <- import("pyteomics")
py_run_string("
import pandas as pd
from pyteomics import mgf
")
source_python("src/snakefiles/functions.py")

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())
PROTON <- 1.007276466622

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  PEAKS_PSM <- fread(snakemake@input[["PEAKS_PSM"]])
  allowed_modifications <- snakemake@input[["allowed_modifications"]]
  
  # directories
  dir_RAW <- snakemake@params[["dir_RAW"]]
  dir_mass_list <- snakemake@params[["dir_mass_list"]]
  dir_RT_calibration <- snakemake@params[["dir_RT_calibration"]]
  
  # RAW file Parser
  ThermoRawFileParser <- "/bin/ThermoRawFileParser/ThermoRawFileParser.exe"
  
  # Resources
  Ncpu <- snakemake@params[["max_cpus"]]
  
} else {
  ### Manual startup
  Experiment_design <- fread("data/Experiment_design.csv") %>% as_tibble()
  PEAKS_PSM <- fread("data/DB search psm.csv")
  # allowed_modifications <- c("Oxidation (M)")
  allowed_modifications <- c("")
  
  # directories
  dir_RAW <- "data/RAW"
  dir_mass_list <- "data/MS_mass_lists/"
  dir_RT_calibration <- "data/RT_calibration/"
  
  # RAW file Parser
  ThermoRawFileParser <- "bin/ThermoRawFileParser.exe"
  
  # Resources
  Ncpu <- parallelly::availableCores() - 1
}

# CPU usage
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl, debug=T)
data.table::setDTthreads(Ncpu)
print(paste0("number of CPUs: ", Ncpu))

# Directories
dir.create(dir_mass_list, recursive = T)
dir.create(dir_RT_calibration, recursive = T)

### ------------------------------------ 1. Convert RAW to MGF --------------------------------------------
# Define paths
input_RAW <- Experiment_design$Filename

output_MGF <- paste0(dir_mass_list, "/", input_RAW, ".mgf")
output_mass_list <- paste0(dir_mass_list, "/", input_RAW, ".txt")
input_RAW <- paste0(dir_RAW, "/", input_RAW, ".raw")

# Prepare ThermoRawFileParser commands
cmds <- paste("mono", ThermoRawFileParser, 
              "-i", input_RAW,
              "-f=0 -l 4",
              "-b", output_MGF)

# Check if files exist - don't repeat RAW -> MGF conversion
if (any(file.exists(output_MGF))) {
  cat(as.character(Sys.time()), " - ", "Found existing MGF files", "\n")
  cmds <- cmds[!file.exists(output_MGF)]
}

if (length(cmds) > 0) {
  cat(as.character(Sys.time()), " - ", "Start ThermoRawFileParser", "\n")
  
  # Convert
  mclapply(cmds, mc.cores = Ncpu, mc.preschedule=T, function(x){
    system(x, intern = T)
  })
} else {
  cat(as.character(Sys.time()), " - ", "All RAW files are converted, nothing to re-run", "\n")
}

### ------------------------------------ 2. Export mass lists --------------------------------------------
lapply(seq_along(output_MGF), function(i){
  cat(as.character(Sys.time()), " - ", "Start reading: ", Experiment_design$Filename[[i]], "\n")
  mgf <- py$generate_mass_rt_list_single_file(output_MGF[[i]]) %>%
    as.data.table() 
  
  # Export mass lists
  cat(as.character(Sys.time()), " - ", "Exporting MW-RT", "\n")
  mass_out <- mgf %>%
    .[, .(precursorMass, retentionTime)] %>%
    setkey(precursorMass, retentionTime) %>%
    unique()
  fwrite(mass_out, file = output_mass_list[[i]], 
         col.names = F,
         append = F, 
         nThread = Ncpu, 
         sep = " ", 
         showProgress = F)
})


### ------------------------------------ 3. Export RT calibration --------------------------------------------
cat(as.character(Sys.time()), " - ", "Exporting PEAKS peptide-RT pairs", "\n")

# Pre-filter the PEAKS PSMs
PEAKS_PSM_filt <- PEAKS_PSM %>%
  # Remove contaminants and decoys
  .[!is.na(Accession)] %>%
  .[!str_detect(Accession, "DECOY")] %>%
  .[,row_ID := .I] %>%
  .[,pep_seq := gsub("[^A-Za-z]", "", Peptide)] # Remove PTM info
    
### Keep only allowed PTMs
if (length(allowed_modifications) > 0) {
  keep_PSMs <- PEAKS_PSM_filt %>%
    .[, .(row_ID, PTM)] %>%
    separate_rows(PTM, sep = ";") %>%
    as.data.table() %>%
    .[PTM %in% c("", allowed_modifications)] 
  
  PEAKS_PSM_filt <- PEAKS_PSM_filt[row_ID %in% keep_PSMs$row_ID,]
}

# Export RT info for calibration
PEAKS_PSM_out <- PEAKS_PSM_filt %>%
  .[, `Source File` := str_remove_all(`Source File`, pattern = ".raw")] %>%
  .[`Source File` %in% Experiment_design$Filename] %>%
  .[, .(`Source File`, pep_seq, RT, Scan)] %>%
  setnames(new = c("Source File", "s", "rt", "scans")) %>%
  # convert PEAKS RT to seconds
  .[, rt := 60 * rt] %>% 
  unique() %>%
  split(by = "Source File", keep.by = FALSE)

# Check if all the files have been searched
if (any(!Experiment_design$Filename %in% names(PEAKS_PSM_out))) {
  missing_files <- Experiment_design$Filename[!Experiment_design$Filename %in% names(PEAKS_PSM_out)]
  
  cat(as.character(Sys.time()), " - ", "Error: missing PEAKS RAW files:", 
      str_c(missing_files, collapse = ";"),
      "\n")
  stop()
}

# Save tables
lapply(seq_along(PEAKS_PSM_out), function(i){
  PEAKS_PSM_out[[i]] %>%
    fwrite(file = paste0(dir_RT_calibration, "/", names(PEAKS_PSM_out)[[i]], ".csv"), 
           append = F, nThread = Ncpu, sep = ",", showProgress = F)
})

### ------------------------------------ Output --------------------------------------------
out <- rbindlist(PEAKS_PSM_out, idcol = "Source File")

if (exists("snakemake")) {
  # fwrite(out, file = unlist(snakemake@output[["PEAKS_PSM_out"]]))
  SPIsnake_log()
  sink()
} else {
  dir.create("results/DB_exhaustive", recursive = T)
  # fwrite(out, file = paste0("results/DB_exhaustive/PEAKS_PSM_out.csv"))
}
