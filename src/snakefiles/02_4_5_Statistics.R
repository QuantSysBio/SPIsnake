### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Compute statistics of SPIsnake databases
#               
# input:        1. Generate_peptides rules are done
#               2. Parameters: which stats to compute
#               3. Experiment_design, Master_table_expanded
# output:       
#               - Statistics of peptide search space: strata sizes and multimapping
#               
# author:       YH, JL, KP

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  # Master_table_expanded
  Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]]) %>% as_tibble()
  
  # Experiment_design
  Experiment_design <- fread(snakemake@input[["Experiment_design"]]) %>% as_tibble()
  
  # CPUs
  Ncpu <- snakemake@params[["max_cpus"]]
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Which stats to make
  Stats_ctrl <- list()
  Stats_ctrl$strata_sizes <- snakemake@params[["strata_sizes"]]
  Stats_ctrl$filtering_sizes <- snakemake@params[["filtering_sizes"]]
  Stats_ctrl$multimapping_within_strata <- snakemake@params[["multimapping_within_strata"]]
  Stats_ctrl$multimapping_between_strata <- snakemake@params[["multimapping_between_strata"]]
  Stats_ctrl$multimapping_filtering_steps <- snakemake@params[["multimapping_filtering_steps"]]
  Stats_ctrl$PTM_stats_masslist <- snakemake@params[["PTM_stats_masslist"]]
  Stats_ctrl$PTM_stats_peptide <- snakemake@params[["PTM_stats_peptide"]]
  
  # Directories
  dir_DB_exhaustive <- snakemake@params[["dir_DB_exhaustive"]]
  dir_DB_PTM_mz <- snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out <- snakemake@params[["dir_DB_out"]]
  
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
  
  # CPUs
  Ncpu <- parallelly::availableCores() - 1
  cl <- parallel::makeForkCluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  data.table::setDTthreads(Ncpu)
  print(paste0("number of CPUs: ", Ncpu))
  retry_times = 3
  
  # Which stats to make
  Stats_ctrl <- list()
  Stats_ctrl$strata_sizes <- TRUE
  Stats_ctrl$filtering_sizes <- TRUE
  Stats_ctrl$multimapping_within_strata <- TRUE
  Stats_ctrl$multimapping_between_strata <- TRUE
  Stats_ctrl$multimapping_filtering_steps <- TRUE
  Stats_ctrl$PTM_stats_masslist <- TRUE
  Stats_ctrl$PTM_stats_peptide <- TRUE  
  
  # Directories
  dir_DB_exhaustive <- "results/DB_exhaustive/"
  dir_DB_PTM_mz <- "results/DB_PTM_mz/"
  dir_DB_out <- "results/DB_out/"
  
  # Wildcard
  filename <- "results/DB_out/.Aggregare_FASTA.done"
  filename_out <- filename
}
SPIsnake_stats <- list()
suppressWarnings(dir.create(dir_DB_out))
suppressWarnings(dir.create(paste0(dir_DB_out, "/Stats/")))

DB_arrow <- open_dataset(paste0(dir_DB_out, "/arrow/"))


### ---------------------------- (2) Strata sizes  --------------------------------------
if (Stats_ctrl$strata_sizes) {
  SPIsnake_stats$strata_sizes <- DB_arrow %>%
    select(peptide, enzyme, MiSl, proteome, index, length) %>%
    group_by(enzyme, MiSl, proteome, index, length) %>%
    summarise(`#unique peptide sequences` = n_distinct(peptide)) %>%
    collect()
  
  strata_sizes <- list()
  strata_sizes$PCP.PSP <- SPIsnake_stats$strata_sizes %>%
    filter(enzyme %in% c("PCP", "PSP")) %>%
    group_by(proteome, enzyme, MiSl, length) %>%
    summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`)) %>%
    mutate(length = as.character(length))
  
  strata_sizes$enzyme <- SPIsnake_stats$strata_sizes %>%
    filter(!enzyme %in% c("PCP", "PSP")) %>%
    group_by(proteome, enzyme, MiSl) %>%
    summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
              length = paste0(min(length), "_", max(length)))
  SPIsnake_stats$strata_sizes <- bind_rows(strata_sizes)
  rm(strata_sizes)
}

### ---------------------------- (3) Strata sizes  ~ filtering steps --------------------------------------
if (Stats_ctrl$filtering_sizes) {
  
  filtering_sizes <- list()
  i <- 1
  for (bg in seq_along(unique(Experiment_design$Biological_group))) {
    Biological_group_bg <- unique(Experiment_design$Biological_group)[[bg]] 
    select_mass_list <- Experiment_design %>%
      filter(Biological_group == Biological_group_bg) %>%
      as_tibble() %>%
      pull(Filename)
    
    ### Molecular weight filter
    pep <- DB_arrow %>%
      dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                    starts_with(paste0("MW.exists:", select_mass_list)), ) %>%
      filter(if_any(starts_with(paste0("MW.exists:", select_mass_list)), ~ .)) %>%
      group_by(enzyme, MiSl, proteome, index, length) %>%
      summarise(`#unique peptide sequences` = n_distinct(peptide),
                Biological_group = Biological_group_bg) %>%
      collect()
    
    pep_MW.PCP.PSP <- pep %>%
      filter(enzyme %in% c("PCP", "PSP")) %>%
      group_by(proteome, enzyme, MiSl, length) %>%
      summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                Biological_group = Biological_group_bg) %>%
      mutate(length = as.character(length),
             filtering_step = "MW.filter")
    
    pep_MW.enzyme <- pep %>%
      filter(!enzyme %in% c("PCP", "PSP")) %>%
      group_by(proteome, enzyme, MiSl) %>%
      summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                Biological_group = Biological_group_bg,
                length = paste0(min(length), "_", max(length)),
                filtering_step = "MW.filter")
    
    filtering_sizes[[i]] <- pep_MW.PCP.PSP
    i <- i + 1    
    filtering_sizes[[i]] <- pep_MW.enzyme
    i <- i + 1
    
    ### Molecular weight + retention time filter
    pep <- DB_arrow %>%
      dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                    starts_with(paste0("MW.RT.exists:", select_mass_list)), ) %>%
      filter(if_any(starts_with(paste0("MW.RT.exists:", select_mass_list)), ~ .)) %>%
      group_by(enzyme, MiSl, proteome, index, length) %>%
      summarise(`#unique peptide sequences` = n_distinct(peptide),
                Biological_group = Biological_group_bg) %>%
      collect()
    
    pep.PCP.PSP <- pep %>%
      filter(enzyme %in% c("PCP", "PSP")) %>%
      group_by(proteome, enzyme, MiSl, length) %>%
      summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                Biological_group = Biological_group_bg) %>%
      mutate(length = as.character(length),
             filtering_step = "MW.RT.filter")
    
    pep.enzyme <- pep %>%
      filter(!enzyme %in% c("PCP", "PSP")) %>%
      group_by(proteome, enzyme, MiSl) %>%
      summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                Biological_group = Biological_group_bg,
                length = paste0(min(length), "_", max(length)),
                filtering_step = "MW.RT.filter")
    
    filtering_sizes[[i]] <- pep.PCP.PSP
    i <- i + 1    
    filtering_sizes[[i]] <- pep.enzyme
    i <- i + 1
    
    ### Molecular weight + retention time + IC50 filter
    alleles <- Experiment_design$`MHC-I_alleles` %>%
      str_c(collapse = "|") %>%
      str_replace_all(pattern = " ", replacement = "")
    
    pep <- DB_arrow %>%
      dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                    starts_with(paste0("Predicted_binder:", select_mass_list)), ) %>%
      filter(if_any(starts_with(paste0("Predicted_binder:", select_mass_list)), ~ .))  %>%
      group_by(enzyme, MiSl, proteome, index, length) %>%
      summarise(`#unique peptide sequences` = n_distinct(peptide),
                Biological_group = Biological_group_bg) %>%
      collect()
    
    pep.PCP.PSP <- pep %>%
      filter(enzyme %in% c("PCP", "PSP")) %>%
      group_by(proteome, enzyme, MiSl, length) %>%
      summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                Biological_group = Biological_group_bg) %>%
      mutate(length = as.character(length),
             filtering_step = "MW.RT.IC50.filter")


    filtering_sizes[[i]] <- pep.PCP.PSP
    i <- i + 1
    rm(pep, pep.PCP.PSP, pep.enzyme)
  }
  SPIsnake_stats$filtering_sizes <- bind_rows(filtering_sizes) 
}

### ---------------------------- (4) Multimapping within strata --------------------------------------
if (Stats_ctrl$multimapping_within_strata) {
  # Unfiltered
  
  if (Stats_ctrl$multimapping_filtering_steps) {
    # MW
    
    # MW.RT
    
    # MW.RT.IC50
  }
}
### ---------------------------- (5) Multimapping between strata --------------------------------------
if (Stats_ctrl$multimapping_between_strata) {
  # Unfiltered
  
  if (Stats_ctrl$multimapping_filtering_steps) {
    # MW
    
    # MW.RT
    
    # MW.RT.IC50
  }
}
### ---------------------------- (6) PTM stats per mass_list --------------------------------------
if (Stats_ctrl$PTM_stats_masslist) {
  
}
### ---------------------------- (7) PTM stats per peptide --------------------------------------
if (Stats_ctrl$PTM_stats_peptide) {
  
}

### ---------------------------- (8) Save outputs --------------------------------------
for (i in seq_along(SPIsnake_stats)) {
  fwrite(SPIsnake_stats[[i]], file = paste0(dir_DB_out, "/Stats/", names(SPIsnake_stats)[[i]], ".csv"), 
         append = F, 
         sep = ",", 
         nThread = Ncpu)
}

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
