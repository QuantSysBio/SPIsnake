### ---------------------------------------------- Aggregate_PTM_mz_RT_matching ----------------------------------------------
# description:  Save SPIsnake main outputs
# 
# output:       1. Arrow: peptide sequences with filter info
#               2. For all filtering steps: FASTA and .CSV
#               
# author:       YH

### ---------------------------- (1) Peptde information: arrow --------------------------------------
### Peptide info
cat("Start saving peptide sequences with filter info", "\n", as.character(Sys.time()), "\n")
if ("chunks" %in% colnames(pep)) {
  pep$chunks <- NULL
}
pep <- setkey(pep, enzyme, index, length)
pep[1:nrow(pep)] %>%  
  as_arrow_table() %>%
  mutate(proteome = chunk_params$proteome) %>%
  mutate(MiSl = chunk_params$MiSl) %>%
  mutate(chunk = chunk_params$chunk) %>%
  relocate(enzyme, MiSl, proteome, index, length, chunk, peptide, MW.unmodified) %>%
  group_by(enzyme, MiSl, proteome, index, length, chunk) %>%
  write_dataset(path = paste0(dir_DB_PTM_mz, "/peptide_seqences/"), 
                existing_data_behavior = "overwrite",
                format = "parquet", 
                compression = "lz4") 
cat("Saved peptide sequences with filter info: Done", "\n", as.character(Sys.time()), "\n")

### ---------------------------- Peptde information: FASTA + CSV --------------------------------------
### ---------------------------- (2) unfiltered --------------------------------------
if (FASTA_outputs_unfiltered){
  cat("Start saving unfiltered peptide sequences", "\n", as.character(Sys.time()), "\n")
  suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/unfiltered")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/unfiltered/chunks/")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/unfiltered")))
  suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/unfiltered/chunks/")))
  
  # Save FASTA
  pep[,.(peptide)] %>%
    .[, fasta := paste0(">", 1:.N, "\n", peptide)] %>%
    .[, .(fasta)] %>%
    vroom_write(file = paste0(dir_DB_out, "FASTA/unfiltered/chunks/", filename, ".fasta"), 
                delim = "\n", col_names = F, append = F, quote = "none", num_threads = Ncpu)
  
  # Save CSV
  vroom_write(pep, append = F, delim = ",", num_threads = Ncpu, 
         file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/unfiltered/chunks/", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
  
  peptides[peptide %chin% pep$peptide,] %>%
    vroom_write(append = F, delim = ",", num_threads = Ncpu, 
           file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/unfiltered/chunks/peptide_mapping_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
  cat("Done saving unfiltered peptide sequences", "\n", as.character(Sys.time()), "\n")
}


### Save filtered peptide sequences
for (bg in seq_along(unique(Experiment_design$Biological_group))) {
  Biological_group_bg <- unique(Experiment_design$Biological_group)[[bg]] 
  select_mass_list <- Experiment_design %>%
    filter(Biological_group == Biological_group_bg) %>%
    as_tibble() %>%
    pull(Filename)
  
  ### ---------------------------- (3) MW filter --------------------------------------
  if (FASTA_outputs_MW_filtered){
    cat("Start saving MW filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW_filtered/chunks")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW_filtered/chunks")))
    
    pep_bg <- pep %>%
      lazy_dt() %>%
      dplyr::select(peptide, index, length, enzyme, MW.unmodified, 
                    starts_with(paste0("MW:", select_mass_list)), 
                    starts_with(paste0("MW.exists:", select_mass_list)), 
                    starts_with(paste0("RT:", select_mass_list)), 
                    starts_with(paste0("MW.RT.exists:", select_mass_list)), 
                    starts_with(paste0("Aff(nM):", select_mass_list)), 
                    starts_with(paste0("Predicted_binder:", select_mass_list, ":"))) %>%
      filter(if_any(starts_with(paste0("MW.exists:", select_mass_list)), ~ .)) %>%
      as.data.table()
    
    if (nrow(pep_bg) > 0){
      # Save FASTA
      pep_bg[,.(peptide)] %>%
        .[, fasta := paste0(">", 1:.N, "\n", peptide)] %>%
        .[, .(fasta)] %>%
        vroom_write(file = paste0(dir_DB_out, "FASTA/MW_filtered/chunks/", Biological_group_bg, "_", filename, ".fasta"), 
                    delim = "\n", col_names = F, append = F, quote = "none", num_threads = Ncpu)
      
      # Save CSV
      vroom_write(pep_bg, append = F, delim = ",", num_threads = Ncpu, 
             file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW_filtered/chunks/", 
                                                     Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
      
      peptides[peptide %chin% pep_bg$peptide,] %>%
        vroom_write(append = F, delim = ",", num_threads = Ncpu, 
               file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW_filtered/chunks/peptide_mapping_", 
                                                       Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))  
    } else {
      cat(paste0("WARNING: NO MW filtered peptides found for Biological_group: ", Biological_group_bg, "!"))
    }
    cat("Done saving MW filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
  }
  rm(pep_bg)
  
  ### ---------------------------- (4) MW.RT filter --------------------------------------
  if (FASTA_outputs_MW_RT_filtered){
    cat("Start saving MW.RT filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW.RT_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW.RT_filtered/chunks")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW.RT_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW.RT_filtered/chunks")))
    
    pep_bg <- pep %>%
      lazy_dt() %>%
      dplyr::select(peptide, index, length, enzyme, MW.unmodified, 
                    starts_with(paste0("MW:", select_mass_list)), 
                    starts_with(paste0("MW.exists:", select_mass_list)), 
                    starts_with(paste0("RT:", select_mass_list)), 
                    starts_with(paste0("MW.RT.exists:", select_mass_list)), 
                    starts_with(paste0("Aff(nM):", select_mass_list)), 
                    starts_with(paste0("Predicted_binder:", select_mass_list, ":"))) %>%
      filter(if_any(starts_with(paste0("MW.RT.exists:", select_mass_list)), ~ .)) %>%
      as.data.table()
    
    if (nrow(pep_bg) > 0){
      pep_bg[,.(peptide)] %>%
        .[, fasta := paste0(">", 1:.N, "\n", peptide)] %>%
        .[, .(fasta)] %>%
        vroom_write(file = paste0(dir_DB_out, "FASTA/MW.RT_filtered/chunks/", Biological_group_bg, "_", filename, ".fasta"), 
                    delim = "\n", col_names = F, append = F, quote = "none", num_threads = Ncpu)
      
      # Save CSV
      vroom_write(pep_bg, append = F, delim = ",", num_threads = Ncpu, 
             file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW.RT_filtered/chunks/", 
                                                     Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
      
      peptides[peptide %chin% pep_bg$peptide,] %>%
        vroom_write(append = F, delim = ",", num_threads = Ncpu, 
               file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW.RT_filtered/chunks/peptide_mapping_", 
                                                       Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))  

    } else {
      cat(paste0("WARNING: NO MW.RT filtered peptides found for Biological_group: ", Biological_group_bg, "!"))
    }
    cat("Done saving MW.RT filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
  }
  
  ### ---------------------------- (5) MW.RT.IC50 filter --------------------------------------
  if (FASTA_outputs_MW_RT_IC50_filtered){
    cat("Start saving MW.RT.IC50 filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW.RT.IC50_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW.RT.IC50_filtered/chunks")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW.RT.IC50_filtered")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW.RT.IC50_filtered/chunks")))
    
    pep_bg <- pep %>%
      lazy_dt() %>%
      dplyr::select(peptide, index, length, enzyme, MW.unmodified, 
                    starts_with(paste0("MW:", select_mass_list)), 
                    starts_with(paste0("MW.exists:", select_mass_list)), 
                    starts_with(paste0("RT:", select_mass_list)), 
                    starts_with(paste0("MW.RT.exists:", select_mass_list)), 
                    starts_with(paste0("Aff(nM):", select_mass_list)), 
                    starts_with(paste0("Predicted_binder:", select_mass_list, ":"))) %>%
      filter(if_any(starts_with(paste0("Predicted_binder:", select_mass_list, ":")), ~ .)) %>%
      as.data.table()
    
    if (nrow(pep_bg) > 0){
      # Save FASTA
      pep_bg[,.(peptide)] %>%
        .[, fasta := paste0(">", 1:.N, "\n", peptide)] %>%
        .[, .(fasta)] %>%
        vroom_write(file = paste0(dir_DB_out, "FASTA/MW.RT.IC50_filtered/chunks/", Biological_group_bg, "_", filename, ".fasta"), 
                    delim = "\n", col_names = F, append = F, quote = "none", num_threads = Ncpu)
      # Save CSV
      vroom_write(pep_bg, append = F, delim = ",", num_threads = Ncpu, 
             file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW.RT.IC50_filtered/chunks/", 
                                                     Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
      
      peptides[peptide %chin% pep_bg$peptide,] %>%
        vroom_write(append = F, delim = ",", num_threads = Ncpu, 
               file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW.RT.IC50_filtered/chunks/peptide_mapping_", 
                                                       Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv"))))) 
    } else {
      cat(paste0("WARNING: NO MW.RT.IC50 filtered peptides found for Biological_group: ", Biological_group_bg, "!"))
    }
    cat("Done saving MW.RT.IC50 filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
  }
  
  ### ---------------------------- (6) MW-PTM filter--------------------------------------
  if (FASTA_outputs_MW_filtered_PTM){
    cat("Start saving MW.PTM filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW_filtered_PTM")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/FASTA/MW_filtered_PTM/chunks")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW_filtered_PTM")))
    suppressWarnings(dir.create(paste0(dir_DB_out, "/CSV/MW_filtered_PTM/chunks")))
    
    # Find peptide chunks
    PTMs <- list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/"), recursive = T)
    PTMs <- PTMs[str_detect(PTMs, paste0("enzyme=", ifelse(chunk_params$enzyme_type == "cis-PSP", "cis-PSP|PSP|PCP", chunk_params$enzyme_type)))]
    PTMs <- PTMs[str_detect(PTMs, paste0("MiSl=", chunk_params$MiSl))]
    PTMs <- PTMs[str_detect(PTMs, paste0("proteome=", chunk_params$proteome))]
    PTMs <- PTMs[str_detect(PTMs, paste0("chunk=", chunk_params$chunk))]
    PTMs <- PTMs[str_detect(PTMs, paste0("mzList=", str_c(select_mass_list, collapse = "|")))]
    
    if (length(PTMs > 0)) {
      pep_bg <- open_dataset(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/")) %>%
        filter(enzyme %in% keep_enzyme) %>%
        filter(MiSl == chunk_params$MiSl) %>%
        filter(proteome == chunk_params$proteome) %>%
        filter(chunk == chunk_params$chunk) %>%
        select(peptide, MW, ids, mzList, enzyme) %>%
        collect() %>%
        setDT()
    } else {
      cat(paste0("WARNING: NO PTMs found for ", Biological_group_bg, "! \n",
                 "Make sure that PTM generation is on for this Biological_group"))
    }
    
    if (nrow(pep_bg) > 0){
      # Save FASTA
      pep_bg[,.(peptide)] %>%
        .[, fasta := paste0(">", 1:.N, "\n", peptide)] %>%
        .[, .(fasta)] %>%
        vroom_write(file = paste0(dir_DB_out, "FASTA/MW_filtered_PTM/chunks/", Biological_group_bg, "_", filename, ".fasta"), 
                    delim = "\n", col_names = F, append = F, quote = "none", num_threads = Ncpu)
      
      # Save CSV
      vroom_write(pep_bg, append = F, delim = ",", num_threads = Ncpu, 
             file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW_filtered_PTM/chunks/", 
                                                     Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv")))))
      
      peptides[peptide %chin% pep_bg$peptide,] %>%
        vroom_write(append = F, delim = ",", num_threads = Ncpu, 
               file = pipe(sprintf("pigz > %s", paste0(dir_DB_out, "CSV/MW_filtered_PTM/chunks/peptide_mapping_", 
                                                       Biological_group_bg, "_", filename, ifelse(compress_CSV, ".csv.gz", ".csv"))))) 
    } else {
      cat(paste0("WARNING: 0-length PTMs found for ", Biological_group_bg, "! \n",
                 "Make sure that PTM generation is on for this Biological_group"))
    }
    cat("Done saving MW.PTM filtered peptide sequences:", Biological_group_bg, "\n", as.character(Sys.time()), "\n")
  }
}

### Clean-up
suppressWarnings(rm(pep_bg, fa, PTMs, bg, Biological_group, Biological_group_bg, select_mass_list))
