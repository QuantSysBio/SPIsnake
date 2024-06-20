### ---------------------------------------------- SPIsnake: stats ----------------------------------------------
# description:  Compute statistics of SPIsnake databases
#               
# input:        1. Generate_peptides rules are done
#               2. Parameters: which stats to compute
#               3. Experiment_design, Master_table_expanded
# output:       
#               - Statistics of peptide search space: strata sizes and multimapping
#               
# author:       YH, JL

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
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
  
  # Directories
  dir_DB_exhaustive <- snakemake@params[["dir_DB_exhaustive"]]
  dir_DB_PTM_mz <- snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out <- snakemake@params[["dir_DB_out"]]
  
  # Wildcard
  filename <- snakemake@output[[1]]
  filename_out <- filename
  
} else {
  ### Manual start
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

### Load arrow database
DB_arrow <- open_dataset(paste0(dir_DB_out, "/arrow/"))

# Define jobs
groups <- DB_arrow$files %>%
  as_tibble() %>%
  mutate(path = value,
         value = str_split_fixed(value, "/arrow/", 2)[,2]) %>%
  mutate(value = str_split_fixed(value, "/part-", 2)[,1]) %>%
  mutate(value = str_remove_all(value, "enzyme=|MiSl=|proteome=|index=|length=")) %>%
  separate(value, into = c("enzyme", "MiSl", "proteome", "index", "length"), sep = "/")  %>%
  mutate(MiSl = as.integer(MiSl),
         length = as.integer(length),
         id = paste0(enzyme, MiSl, proteome, length, index))

print(glimpse(groups))

### ---------------------------- (2) Strata sizes  --------------------------------------
if (Stats_ctrl$strata_sizes) {
  cat(as.character(Sys.time()), " - ", "Estimate strata sizes: Starting", "\n")
  
  SPIsnake_stats$strata_sizes <- DB_arrow %>%
    group_by(index, length, proteome, enzyme, MiSl) %>%
    summarise(`#unique peptide sequences` = n()) %>%
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
  cat(as.character(Sys.time()), " - ", "Estimate strata sizes: Done", "\n")
}

### ---------------------------- (3) Strata sizes  ~ filtering steps --------------------------------------
if (Stats_ctrl$filtering_sizes) {
  cat(as.character(Sys.time()), " - ", "Estimate filtered strata sizes: Starting", "\n")
  
  filtering_sizes <- list()
  i <- 1
  
  ### By mass_list - lower RAM usage
  for (bg in seq_along(unique(Experiment_design$Filename))) {
    Biological_group_bg <- unique(Experiment_design$Filename)[[bg]] 
    select_mass_list <- Experiment_design %>%
      filter(Experiment_design$Filename == Biological_group_bg) %>%
      as_tibble() %>%
      pull(Filename)
    print(Biological_group_bg)
    print(select_mass_list)
    
    ### Molecular weight filter
    if (any(str_starts(names(open_dataset(DB_arrow$files[[1]])), "MW.exists:"))) {
      pep <- DB_arrow %>%
        dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                      starts_with(paste0("MW.exists:", select_mass_list))) %>%
        group_by(enzyme, MiSl, proteome, index, length) %>%
        filter(if_any(starts_with(paste0("MW.exists:", select_mass_list)), ~ .)) %>%
        summarise(`#unique peptide sequences` = n(),
                  Biological_group = Biological_group_bg) %>%
        collect()
      
      pep_MW.PCP.PSP <- pep %>%
        filter(enzyme %in% c("PCP", "PSP")) %>%
        group_by(proteome, enzyme, MiSl, length) %>%
        summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                  Biological_group = Biological_group_bg) %>%
        mutate(length = as.character(length),
               filtering_step = "mz_matched")
      
      pep_MW.enzyme <- pep %>%
        filter(!enzyme %in% c("PCP", "PSP")) %>%
        group_by(proteome, enzyme, MiSl) %>%
        summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                  Biological_group = Biological_group_bg,
                  length = paste0(min(length), "_", max(length)),
                  filtering_step = "mz_matched")
      
      filtering_sizes[[i]] <- pep_MW.PCP.PSP
      i <- i + 1    
      filtering_sizes[[i]] <- pep_MW.enzyme
      i <- i + 1
      cat(as.character(Sys.time()), " - ", "Molecular weight: Done", "\n")
    }
    
    ### Molecular weight + retention time filter
    if (any(str_starts(names(open_dataset(DB_arrow$files[[1]])), "MW.RT.exists:"))) {
      pep <- DB_arrow %>%
        dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                      starts_with(paste0("MW.RT.exists:", select_mass_list)), ) %>%
        group_by(enzyme, MiSl, proteome, index, length) %>%
        filter(if_any(starts_with(paste0("MW.RT.exists:", select_mass_list)), ~ .)) %>%
        summarise(`#unique peptide sequences` = n(),
                  Biological_group = Biological_group_bg) %>%
        collect()
      
      pep.PCP.PSP <- pep %>%
        filter(enzyme %in% c("PCP", "PSP")) %>%
        group_by(proteome, enzyme, MiSl, length) %>%
        summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                  Biological_group = Biological_group_bg) %>%
        mutate(length = as.character(length),
               filtering_step = "mz_RT_matched")
      
      pep.enzyme <- pep %>%
        filter(!enzyme %in% c("PCP", "PSP")) %>%
        group_by(proteome, enzyme, MiSl) %>%
        summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                  Biological_group = Biological_group_bg,
                  length = paste0(min(length), "_", max(length)),
                  filtering_step = "mz_RT_matched")
      
      filtering_sizes[[i]] <- pep.PCP.PSP
      i <- i + 1    
      filtering_sizes[[i]] <- pep.enzyme
      i <- i + 1
      cat(as.character(Sys.time()), " - ", "Molecular weight + retention time: Done", "\n")
    }
    
    ### Molecular weight + retention time + IC50 filter
    if (any(str_starts(names(open_dataset(DB_arrow$files[[1]])), "Predicted_binder:"))) {
      alleles <- Experiment_design$`MHC-I_alleles` %>%
        str_c(collapse = "|") %>%
        str_replace_all(pattern = " ", replacement = "")
      
      pep <- DB_arrow %>%
        dplyr::select(peptide, enzyme, MiSl, proteome, index, length,
                      starts_with(paste0("Predicted_binder:", select_mass_list))) %>%
        group_by(enzyme, MiSl, proteome, index, length) %>%
        filter(if_any(starts_with(paste0("Predicted_binder:", select_mass_list)), ~ .))  %>%
        summarise(`#unique peptide sequences` = n(),
                  Biological_group = Biological_group_bg) %>%
        collect()
      
      pep.PCP.PSP <- pep %>%
        filter(enzyme %in% c("PCP", "PSP")) %>%
        group_by(proteome, enzyme, MiSl, length) %>%
        summarise(`#unique peptide sequences` = sum(`#unique peptide sequences`),
                  Biological_group = Biological_group_bg) %>%
        mutate(length = as.character(length),
               filtering_step = "IC50_filtered")
      
      
      filtering_sizes[[i]] <- pep.PCP.PSP
      i <- i + 1
      cat(as.character(Sys.time()), " - ", "Molecular weight + retention time + IC50 filter: Done", "\n")
      rm(pep, pep.PCP.PSP, pep.enzyme)
    }
  }
  SPIsnake_stats$filtering_sizes <- bind_rows(filtering_sizes) 
  SPIsnake_stats$Summary_stats <- SPIsnake_stats$filtering_sizes %>%
    pivot_wider(id_cols = c("proteome", "enzyme", "MiSl", "length", "Biological_group"), 
                names_from = "filtering_step", 
                values_from = "#unique peptide sequences") %>%
    left_join(SPIsnake_stats$strata_sizes) %>%
    rename(Proteome = proteome,
           enzyme_type = enzyme,
           length = length,
           unique = `#unique peptide sequences`) %>%
    relocate(Biological_group, enzyme_type, MiSl, Proteome, length, unique, mz_matched, mz_RT_matched)
  
  cat(as.character(Sys.time()), " - ", "Estimate filtered strata sizes: Done", "\n")
}

### ---------------------------- Save outputs --------------------------------------
### tables
for (i in seq_along(SPIsnake_stats)) {
  fwrite(SPIsnake_stats[[i]], file = paste0(dir_DB_out, "/Stats/", names(SPIsnake_stats)[[i]], ".csv"), 
         append = F, 
         sep = ",", 
         nThread = Ncpu)
}

### plots
### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
