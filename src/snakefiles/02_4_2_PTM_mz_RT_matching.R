### ---------------------------------------------- Aggregate_PTM_mz_RT_matching ----------------------------------------------
# description:  Find a unique set of peptide sequences, compute molecular weight (MW) and do m/z and RT matching with all input mass_lists. 
#               If PTMs are required, generate them too. 
#               
# input:        1. Peptide sequences
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .parquet
#               - Unique peptides after m/z and RT matching per experiment .parquet
#               
# author:       YH, JL, KP

if (method %in% c("achrom", "AutoRT")) {
  use_condaenv("R_env_reticulate")
  pyteomics <- import("pyteomics")
  
  py_run_string("
import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None
")
  
  py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
") 
}
achrom_predict_RT <- function(x, 
                              RCs = RCs){
  if (length(x) == 1) {
    out <- unique(py_calls$achrom_calculate_RT(rep(x, 2), RCs, r_to_py(FALSE)))
  } else {
    out <- py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE))
  }
  return(out)
}

# Define IC50 aggregation
IC50_aggregation_table <- Experiment_design %>%
  dplyr::select(Filename, Biological_group, `MHC-I_alleles`, Affinity_threshold) %>%
  filter(!(is.na(`MHC-I_alleles`) | is.na(Affinity_threshold))) %>%
  tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]", ) %>%
  mutate(`MHC-I_alleles` = str_squish(`MHC-I_alleles`),
         Affinity_threshold = str_squish(Affinity_threshold)) %>%
  unique() 
IC50_aggregation_table

### ---------------------------- (2) Compute MW & PTMS--------------------------------------
### Define fixed_mods ~ mass_lists
mods_fixed <- unique(na.omit(Experiment_design$modifications_fixed))
mods_fixed <- mods_fixed[!mods_fixed == ""]

Experiment_design_fixed_mods <- Experiment_design %>%
  mutate(modifications_fixed = ifelse(modifications_fixed == "", NA, modifications_fixed)) %>%
  select(Filename, modifications_fixed) %>%
  rename(mzList_fixed_mods = Filename)

fixed_mods <- lapply(Experiment_design_fixed_mods$modifications_fixed, function(x){
  if (is.na(x)) {
    tibble(Id = "none",
           Site = "none",
           Position = "none",
           MonoMass = 0)
  } else {
    fread(paste0("data/modifications_fixed/", x, ".csv")) %>% as_tibble()
  }
})
names(fixed_mods) <- Experiment_design_fixed_mods$mzList_fixed_mods

### Define variable_mods ~ proteomes
PTM_list <- left_join(params, Master_table_expanded) %>%
  # mutate(PTMs = "Minimal") %>% ### TODO remove "mutate(PTMs = "Minimal")"
  pull(PTMs) %>%
  str_split(pattern = "\\|") %>%
  unlist() %>%
  str_squish()
PTM_list <- PTM_list[!(PTM_list == "")] %>% unique() %>% na.omit() %>% as.character()

### Whether to generate PTMs for spliced peptides
if (chunk_params$enzyme_type %in% c("PSP", "cis-PSP")) {
  # PCP peptides
  pep_PCP <- peptides[enzyme == "PCP",.(index, peptide, enzyme)] %>%
    unique()
} 
if (chunk_params$enzyme_type %in% c("PSP", "cis-PSP")) {
  if (generate_spliced_PTMs == FALSE) {
    ignore_generate_spliced_PTMs = !(enzyme_type == "PSP" | enzyme_type == "cis-PSP")
    pep <- pep_PCP
    
    print("PTMs will NOT be generated for cis-PSP")
  } else {
    ignore_generate_spliced_PTMs = TRUE
    print("PTMs will be generated for cis-PSP")
    
    # All unique peptides
    pep <- peptides[, .(index, peptide, length, enzyme)]  %>%
      unique() 
  }
} else {
  # All unique peptides
  pep <- peptides[, .(index, peptide, length, enzyme)]  %>%
    unique() 
}

### Define chunking strategy
n_chunks = ceiling((nrow(pep) / PTM_chunk))
cat(as.character(Sys.time()), " - ", "Creating: ", n_chunks, " chunks", "\n")

### Compute MW chunk-wise
pep <- pep %>%
  .[, chunks := rep(1:n_chunks, each=PTM_chunk, length.out = nrow(pep))] %>%
  .[, MW := computeMZ_biostrings(peptide), by = chunks]
cat(as.character(Sys.time()), " - ", "Computed MW for: ", n_chunks, " chunks", "\n")

### PTMs: fixed and/or variable
if (length(PTM_list) > 0 & length(mods_fixed) > 0) {
  PTM_mode = "Fixed_&_variable"
} else if (length(PTM_list) > 0 & length(mods_fixed) == 0) {
  PTM_mode = "Variable_only"
} else if (length(PTM_list) == 0 & length(mods_fixed) > 0) {
  PTM_mode = "Fixed_only"
} else if (length(PTM_list) == 0 & length(mods_fixed) == 0) {
  PTM_mode = "No_mods"
}
print(paste0("PTM mode: ", PTM_mode))

if (PTM_mode != "No_mods") {
  Ncpu_PTM <- ceiling(Ncpu / 2)
  for (PTM in PTM_list) {
    ### Load PTM table
    if (PTM_mode == "Fixed_only") {
      mods <- tibble(Id = "none",
                     Site = "none",
                     Position = "none",
                     MonoMass = 0)
    } else {
      mods <- fread(paste0("data/modifications/", PTM, ".csv"), nThread = Ncpu) %>% as_tibble()
    }
    
    # Prepare peptides
    for (input_i in 1:n_chunks) {
      print(paste("Generating PTMs for:", PTM, input_i))
      
      ### Additional split to reduce RAM usage
      PTMcombinations <- pep[chunks == input_i, .(peptide, MW)]
      if (nrow(PTMcombinations) == 1) {
        PTMcombinations <- rbind(PTMcombinations, PTMcombinations)
      }
      
      if (length(mods_fixed) > 0) {
        
        ### (1) Add fixed mods per mass_list
        fm_out <- list()
        fm_seq <- list()
        for (fm in names(fixed_mods)) {
          fm_PTMs <- select(fixed_mods[[fm]], Site, Position)
          
          fa = filter(fm_PTMs, fm_PTMs$Position=="Anywhere") %>% pull(Site)
          fc = filter(fm_PTMs, Site=="C-term" & fm_PTMs$Position=="Any C-term") %>% pull(Site)
          fn = filter(fm_PTMs, Site=="N-term" & fm_PTMs$Position=="Any N-term") %>% pull(Site)
          
          PTMcombinations_fm <- list()
          if (length(fa) > 0) {
            PTMcombinations_fm[["fa"]] <- PTMcombinations %>%
              filter(str_detect(peptide, str_c(fa, sep = "|"))) %>%
              as.data.table()
          } 
          if (length(fc) > 0) {
            PTMcombinations_fm[["fc"]] <- PTMcombinations %>%
              filter(str_ends(peptide, str_c(fc, sep = "|"))) %>%
              as.data.table()
          } 
          if (length(fn) > 0) {
            PTMcombinations_fm[["fa"]] <- PTMcombinations %>%
              filter(str_starts(peptide, str_c(fn, sep = "|"))) %>%
              as.data.table()
          } 
          PTMcombinations_fm <- rbindlist(PTMcombinations_fm)
          fm_seq[[fm]] <- PTMcombinations_fm
          if (nrow(PTMcombinations_fm) > 0) {
            PTMcombinations_fm <- PTMcombinations_fm[, PTMcombinations := parallel::mcmapply(getPTMcombinations_fixed_vec, 
                                                                                             peptide, MW, 
                                                                                             list(fixed_mods[[fm]]), 
                                                                                             SIMPLIFY = T, 
                                                                                             mc.cores = Ncpu_PTM, 
                                                                                             mc.preschedule = T, 
                                                                                             USE.NAMES = F, 
                                                                                             mc.cleanup = T)]
            PTMcombinations_fm <- rbindlist(PTMcombinations_fm$PTMcombinations, use.names = F) 
            
            ### Generate variable mods per mass list in peptides with fixed mods
            mods_fm <- anti_join(mods, fm_PTMs)
            if (PTM_mode == "Fixed_&_variable") {
              PTMcombinations_fm2 <- PTMcombinations_fm[, PTMcombinations := parallel::mcmapply(getPTMcombinations_fast_vec, 
                                                                                                peptide, MW, max_variable_PTM,
                                                                                                list(mods_fm), 
                                                                                                SIMPLIFY = T, 
                                                                                                mc.cores = Ncpu_PTM, 
                                                                                                mc.preschedule = T, 
                                                                                                USE.NAMES = F, 
                                                                                                mc.cleanup = T)]
              PTMcombinations_fm2 <- PTMcombinations_fm2$PTMcombinations
              names(PTMcombinations_fm2) <- PTMcombinations_fm$ids
              PTMcombinations_fm2 <- rbindlist(PTMcombinations_fm2, use.names = T, idcol = "fixed_mods")  
              PTMcombinations_fm2 <- PTMcombinations_fm2[,ids:=paste(fixed_mods, ids, sep = ";")] %>%
                .[,fixed_mods:=NULL] 
              
              fm_out[[fm]] <- PTMcombinations_fm2
              rm(PTMcombinations_fm, PTMcombinations_fm2)
            } else {
              fm_out[[fm]] <- PTMcombinations_fm
              rm(PTMcombinations_fm)
            }
          }
        } # end fixed mod ~ mass_list
        
        ### (2) Find all peptides without fixed mods
        fm_seq <- fm_seq[unlist(lapply(fm_seq, nrow)) > 0]
        if (length(fm_seq) > 0) {
          vm_out <- PTMcombinations[!peptide %chin% rbindlist(fm_seq)$peptide,]
          if (nrow(vm_out) > 0 ) {
            vm_out <- vm_out[, PTMcombinations := parallel::mcmapply(getPTMcombinations_fast_vec, 
                                                                     peptide, MW, max_variable_PTM, 
                                                                     list(mods), 
                                                                     SIMPLIFY = T, 
                                                                     mc.cores = Ncpu_PTM, 
                                                                     mc.preschedule = T, 
                                                                     USE.NAMES = F, 
                                                                     mc.cleanup = T)]
            vm_seq <- vm_out[,.(peptide, MW)]
          } else {
            vm_seq <- vm_out[,.(peptide, MW)]
          }
        }
        
        ### (3) Find all peptides without fixed mods not common between mass lists
        ### Generate their variable PTMs
        uv_out <- list()
        uv_seq <- list()
        for (fm in names(fixed_mods)) {
          
          if (exists("vm_seq")) {
            uv_out[[fm]] <- PTMcombinations[!peptide %chin% vm_seq$peptide,]
          } else {
            uv_out[[fm]] <- PTMcombinations
          }
          uv_out[[fm]] <- uv_out[[fm]][!peptide %chin% fm_seq[[fm]]$peptide,]
          
          if(nrow(uv_out[[fm]]) > 0){
            uv_out[[fm]] <- uv_out[[fm]][, PTMcombinations := parallel::mcmapply(getPTMcombinations_fast_vec, 
                                                                                 peptide, MW, max_variable_PTM, 
                                                                                 list(mods), 
                                                                                 SIMPLIFY = T, 
                                                                                 mc.cores = Ncpu_PTM, 
                                                                                 mc.preschedule = T, 
                                                                                 USE.NAMES = F, 
                                                                                 mc.cleanup = T)]
            uv_seq[[fm]] <- uv_out[[fm]][,.(peptide, MW)]
          } else {
            uv_out[[fm]] <- NULL
          }
        }
        ### (4) Find all peptides without fixed mods not common between mass lists
        # 4.1 fixed ~ mass list
        PTMcombinations_list <- list()
        PTMcombinations_list$fm <- rbindlist(fm_out, idcol = "mass_list_tag")
        
        # 4.2 unique ~ mass list
        uv <- rbindlist(uv_out, idcol = "mass_list_tag")$PTMcombinations
        names(uv) <- rbindlist(uv_out, idcol = "mass_list_tag")$mass_list_tag
        PTMcombinations_list$uv <- rbindlist(uv, idcol = "mass_list_tag")
        
        # 4.3 all variable mods common across mass lists
        if (exists("vm_out")) {
          vm <- vm_out$PTMcombinations
          vm <- rbindlist(vm)
          vm[,mass_list_tag:="common_variable_mods"]
          PTMcombinations_list$vm <- vm
        }
        
        PTMcombinations_list <- PTMcombinations_list[unlist(lapply(PTMcombinations_list, ncol)) == 4]
        PTMcombinations <- rbindlist(PTMcombinations_list, use.names = T) 
        PTMcombinations <- PTMcombinations %>%
          .[!is.na(MW),] %>%
          .[, id := rowid(peptide, ids)] %>%
          .[, c("MW_Min", "MW_Max") := MW]
        
        #### end fixed mods 
      } else { 
        # No fixed, only variable mods
        PTMcombinations <- PTMcombinations[, PTMcombinations := parallel::mcmapply(getPTMcombinations_fast_vec, 
                                                                                   peptide, MW, max_variable_PTM, 
                                                                                   list(mods), 
                                                                                   SIMPLIFY = T, 
                                                                                   mc.cores = Ncpu_PTM, 
                                                                                   mc.preschedule = T, 
                                                                                   USE.NAMES = F, 
                                                                                   mc.cleanup = T)]
        print("got PTM combinations")
        
        PTMcombinations <- rbindlist(PTMcombinations$PTMcombinations, use.names = F) 
        PTMcombinations <- PTMcombinations %>%
          .[!is.na(ids)] %>%
          .[, id := rowid(peptide, ids)] %>%
          .[, c("MW_Min", "MW_Max") := MW] %>%
          .[,mass_list_tag:="common_variable_mods"]
        print("rbind combinations")
      }
      
      # If there have been PTMs generated
      if (!(nrow(PTMcombinations) == 1 & anyNA(PTMcombinations$ids))) {
        PTM_MW_out <- foverlaps(x = PTMcombinations, y = mzList, 
                                by.x=c("MW_Min","MW_Max"),
                                by.y=c("MW_Min", "MW_Max"), nomatch=0L) %>%
          .[, c("i.MW_Min", "i.MW_Max", "MW_Min", "MW_Max"):=NULL] %>%
          unique() %>%
          .[, id := NULL] %>%
          # Keep combinations where mzList matches mass_list_tag or common variable mods
          .[(mzList == mass_list_tag) | (mass_list_tag == "common_variable_mods"),] 
        print("MW filter: Done")
        
        # Save modified peptides
        if (chunk_params$enzyme_type %in% c("PSP", "cis-PSP")) {
          if (generate_spliced_PTMs == TRUE) {
            PTM_MW_out <- PTM_MW_out[, enzyme := fifelse(peptide %chin% pep_PCP$peptide, "PCP", "PSP")]
          } else {
            PTM_MW_out$enzyme <- "PCP"
          }
        } else {
          PTM_MW_out$enzyme <- chunk_params$enzyme_type
        }
        
        ### Save outputs
        PTM_MW_out[1:nrow(PTM_MW_out)] %>%
          as_arrow_table() %>%
          mutate(Variable_PTMs = PTM) %>%
          mutate(proteome = chunk_params$proteome) %>%
          mutate(MiSl = chunk_params$MiSl) %>%
          mutate(chunk = chunk_params$chunk) %>%
          mutate(input_i = input_i) %>%
          group_by(enzyme, MiSl, proteome, chunk, Variable_PTMs, mzList, input_i) %>%
          write_dataset(path = paste0(dir_DB_PTM_mz, "/unique_peptides_mz_matched_PTM/"), 
                        existing_data_behavior = "overwrite",
                        format = "parquet", 
                        compression = "lz4") 
      } 
      rm(PTMcombinations, PTM_MW_out)
    } # end input_i
  } # end only variable PTM
} # end PTM generation

### clean-up
suppressWarnings(rm(j, pep_PCP, mzList, PTMcombinations_list, PTMcombinations_fm, PTMs, PTM_stats_all,
                    uv, uv_out, uv_seq, vm, vm_out, vm_seq, fm_out, fm_seq, fm_PTMs, fa, fc, fm, fn, input_i))

if (chunk_params$enzyme_type %in% c("PSP", "cis-PSP")) {
  if (generate_spliced_PTMs == FALSE) {
    
    # Keep working with all unique peptides including cis-PSP
    pep <- peptides[, .(index, peptide, length, enzyme)] %>%
      unique() 
    n_chunks = ceiling((nrow(pep) / PTM_chunk))
    
    ### Compute MW chunk-wise
    pep <- pep %>%
      .[, chunks := rep(1:n_chunks, each=2*PTM_chunk, length.out = nrow(pep))] %>%
      .[, MW := computeMZ_biostrings(peptide), by = chunks]
    
    cat(as.character(Sys.time()), " - ", "Computing MW for: PCP and cis-PCP: Done", "\n")
  } 
} 

### ---------------------------- (3) MW & RT matching --------------------------------------
list_IC50 <- list()

# Copy peptides into a new data.table
pep[, MW.unmodified := MW]
pep <- pep[, .(index, peptide, length, enzyme, chunks, MW.unmodified, MW)]

for (j in 1:nrow(MS_mass_lists)) {
  print(MS_mass_lists$mass_list[j])
  MS_mass_list <- MS_mass_lists$mass_list[j]
  MW.exists.j <- paste0("MW.exists:", MS_mass_list)
  MW.RT.exists.j <- paste0("MW.RT.exists:", MS_mass_list)
  
  # Filter tolerances
  tolerance <- as.numeric(Experiment_design$Precursor_mass_tolerance_ppm[Experiment_design$Filename == MS_mass_list])
  RT_tolerance <- as.numeric(RT_Performance_df$mean_value[RT_Performance_df$dataset == MS_mass_list & RT_Performance_df$metric == "MAE"])
  
  # Calibration input
  mzList <- MS_mass_lists_data[[j]]
  mzList <- setkey(mzList, MW_Min, MW_Max, RT_Min, RT_Max)
  
  # If there are fixed mods - re-compute MW for these peptides
  mods_fixed_j <- fixed_mods[[MS_mass_lists$mass_list[j]]]
  if (!mods_fixed_j$Id[[1]] == "none") {
    cat(as.character(Sys.time()), " - ", "Updating masses for fixed PTMs: Starting \n")
    pep[, MW := MW.unmodified]
    
    if (all(mods_fixed_j$Position == "Anywhere")) {
      cat(as.character(Sys.time()), " - ", "All fixed mods Position == Anywhere, using Biostrings for MW computation \n")
      masses_j <- monoisotopic_masses
      for (m in seq_along(mods_fixed_j$Site)) {
        masses_j$mass[masses_j$AA == mods_fixed_j$Site[[m]]] <- masses_j$mass[masses_j$AA == mods_fixed_j$Site[[m]]] + mods_fixed_j$MonoMass[[m]] 
      }
      pep <- pep[str_detect(peptide, str_c(mods_fixed_j$Site, collapse = "|")), 
                 MW := computeMZ_biostrings(peptide, AAs = masses_j$AA, masses = masses_j$mass), by = chunks]
      rm(masses_j)
    } else {
      cat(as.character(Sys.time()), " - ", "Using position-aware fixed PTM MW computation \n")
      pep <- pep[str_detect(peptide, str_c(mods_fixed_j$Site, collapse = "|")),
                 MW := parallel::mcmapply(getPTMcombinations_fixed_mass_vec,
                                          peptide, MW.unmodified,
                                          list(mods_fixed_j),
                                          SIMPLIFY = T,
                                          mc.cores = Ncpu,
                                          mc.preschedule = T,
                                          USE.NAMES = F,
                                          mc.cleanup = T), by = chunks]
    }
    cat(as.character(Sys.time()), " - ", "Updating masses for fixed PTMs: Done \n")
  } else {
    pep[, MW := MW.unmodified]
    cat(as.character(Sys.time()), " - ", "No fixed PTMs: masses will not be updated \n")
  }
  
  # MW filter block-wise
  cat(as.character(Sys.time()), " - ", "MW filter: Starting \n")
  tmp_1D <- pep[MW %inrange% mzList[,c("MW_Min", "MW_Max")], .(peptide), by = chunks]
  pep[, get("MW.exists.j") := fifelse(peptide %chin% tmp_1D$peptide, TRUE, FALSE)]
  cat(as.character(Sys.time()), " - ", "MW filter: Done \n")
  
  ### ----------------------------- Predict RT for m/z matched peptides -----------------------------
  if (method == "AutoRT") {
    cat(as.character(Sys.time()), " - ", "RT prediction: Starting AutoRT \n")
    
    # Save peptides for AutoRT prediction
    pep <- pep[, predict_RT := pep[[get("MW.exists.j")]]]
    
    pep[pep[[get("MW.exists.j")]],] %>%
      rename(x=peptide) %>%
      select(x) %>%
      as.data.table() %>%
      fwrite(file = paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", filename, "_", MS_mass_list, ".tsv"), 
             nThread = Ncpu, append = F, sep = "\t")
    
    # AutoRT predict with pre-trained model
    system(command = paste("python bin/AutoRT/autort.py predict --test",
                           paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", filename, "_", MS_mass_list, ".tsv"),
                           "-s", paste0("results/RT_prediction/RT_models/", MS_mass_list, "/model.json"),
                           "-o", paste0("results/RT_prediction/peptide_RT/", filename, "_", MS_mass_list)),
           intern = T)
    
    # Read in AutoRT predictions
    RT_pred <- fread(paste0("results/RT_prediction/peptide_RT/", filename, "_", MS_mass_list, "/test.tsv"), 
                     nThread = Ncpu)
    pep[pep[[get("MW.exists.j")]], RT_pred := RT_pred$y_pred]
    
    ### Clean-up
    rm(RT_pred)
    file.remove(paste0("results/DB_PTM_mz/unique_peptides_mz_matched/", filename, "_", MS_mass_list, ".tsv"))
    cat(as.character(Sys.time()), " - ", "RT prediction: done AutoRT \n")
    
  } else if (method == "achrom") {
    cat(as.character(Sys.time()), " - ", "RT prediction: Loading parameters \n")
    # Load parameters
    RCs <- readRDS(paste0("results/RT_prediction/RT_models/", MS_mass_list, "/achrom_RCs.rds"))
    
    # Use Leu values for Ile in case no Ile was provided in calibration
    if (is.null(RCs$aa$I)) {
      RCs$aa$I <- RCs$aa$L
    }
    
    # Exclusion pattern: peptides containing these letters will be omitted
    exclusion_pattern <- AA[!AA %in% names(RCs$aa)] %>% str_c(collapse = "|")
    if (!exclusion_pattern == "") {
      # pep <- pep[, predict_RT := pep[[get("MW.exists.j")]] &
      #              str_detect(peptide, pattern = exclusion_pattern, negate = T)]
      pep <- pep[, achrom_exclude := str_detect(peptide, pattern = exclusion_pattern, negate = T), by = chunks]
      pep <- pep[, predict_RT := pep[[get("MW.exists.j")]] & achrom_exclude]
      pep[, achrom_exclude := NULL]
    } else {
      pep <- pep[, predict_RT := pep[[get("MW.exists.j")]]]
    }
    
    ### Predict RT
    not_empty_MS_mass_list <- length(pep$predict_RT) > 0
    if (not_empty_MS_mass_list) {
      cat(as.character(Sys.time()), " - ", "RT prediction: Starting achrom \n")
      pep <- pep[predict_RT==TRUE, RT_pred := split(pep[predict_RT==TRUE,]$peptide, factor(ceiling(seq_along(pep[predict_RT==TRUE,]$peptide)/ceiling(length(pep[predict_RT==TRUE,]$peptide)/Ncpu)))) %>%
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
  cat(as.character(Sys.time()), " - ", "RT prediction: Done \n")
  
  ### 2D filter: MW & RT
  if (not_empty_MS_mass_list) {
    cat(as.character(Sys.time()), " - ", "2D MW.RT filter: Starting \n")
    pep <- setkey(pep, MW, RT_pred)
    tmp_2D <- mzList[pep[predict_RT==TRUE, .(peptide, MW, RT_pred)], 
                     on=.(MW_Min <= MW, 
                          MW_Max >= MW,
                          RT_Min <= RT_pred, 
                          RT_Max >= RT_pred), 
                     nomatch=0, .(peptide)] %>%
      unique()
    
    ### 2D filter: add info to the main table
    pep[, get("MW.RT.exists.j") := FALSE]
    pep[predict_RT==TRUE, get("MW.RT.exists.j") := fifelse(peptide %chin% tmp_2D$peptide, TRUE, FALSE)]
    cat(as.character(Sys.time()), " - ", "2D MW.RT filter: Done \n")
  } else {
    tmp_2D <- data.table()
  }
  setnames(pep, old = "MW", new = paste0("MW:", MS_mass_list))
  setnames(pep, old = "RT_pred", new = paste0("RT:", MS_mass_list))
  if ("predict_RT" %in% colnames(pep)) {
    pep[, predict_RT := NULL]
  }
  
  ### ---------------------------- Define IC50 --------------------------------------
  IC50_aggregation_table_j <- IC50_aggregation_table %>%
    dplyr::filter(Filename == MS_mass_lists$mass_list[j]) %>%
    as_tibble()
  
  if (not_empty_MS_mass_list & 
      nrow(tmp_2D) > 0 & 
      nrow(IC50_aggregation_table_j) > 0 &
      chunk_params$enzyme_type %in% c("PCP", "PSP", "cis-PSP")) {
    
    # Which peptides to predict IC50 for
    Aff_name <- unique(paste0("Aff(nM):", IC50_aggregation_table_j$`MHC-I_alleles`)) 
    pep[, get("Aff_name") := fifelse(.SD == T, T, F), .SDcols = get("MW.RT.exists.j")]
  }
  rm(tmp_2D, not_empty_MS_mass_list)
} # End 2D filter ~ datasets


### ---------------------------- (4) IC50 prediction --------------------------------------
alleles <- IC50_aggregation_table %>% 
  pull(`MHC-I_alleles`) %>%
  unique() %>% 
  na.omit()
Affinity_cols <- paste0("Aff(nM):", alleles)
cat(as.character(Sys.time()), " - ", "Predicting MHC-I affinity for: ", str_c(alleles, collapse = "; "), "\n")

# Are there any sequences for prediction?
if (length(alleles) > 0) {
  peptide_exports <- pep %>%
    lazy_dt() %>%
    select(peptide, Affinity_cols) %>%
    filter(if_any(Affinity_cols, ~ .)) %>%
    as.data.table()
  peptide_exports
} else {
  peptide_exports <- data.table()
}

if (nrow(peptide_exports) > 0 & 
    nrow(IC50_aggregation_table) > 0 &
    chunk_params$enzyme_type %in% c("PCP", "PSP", "cis-PSP")) {
  
  # Export peptide sequences for MHC-I affinity prediction
  peptide_exports_files <- lapply(1:length(Affinity_cols), function(i){
    print(i)
    print(Affinity_cols[[i]])
    out <- peptide_exports %>%
      .[(get(Affinity_cols[[i]])), .(peptide)] %>%
      lazy_dt() %>%
      mutate(n_netMHCpan_chunks = max(ceiling(n() / netMHCpan_chunk), Ncpu)) %>%
      mutate(netMHCpan_split = rep(1:n_netMHCpan_chunks, each=ceiling(n()/n_netMHCpan_chunks), length.out = n())) %>%
      mutate(n_netMHCpan_chunks = NULL) %>%
      as.data.table()
    out <- split(out, by = "netMHCpan_split", drop = T, keep.by = T) 
    
    out <- out %>%
      lapply(FUN = function(x){
        netMHCpan_pep <- paste0(dir_DB_PTM_mz, "/netMHCpan_input/", 
                                filename, "_", 
                                alleles[[i]], "_ch_", x$netMHCpan_split[[1]] ,".tsv")
        
        x[, .(peptide)] %>% 
          fwrite(file = netMHCpan_pep,
                 sep = "\t", nThread = Ncpu, append = FALSE, col.names = FALSE)          
        return(netMHCpan_pep)
      })
    out <- data.table(input = unlist(out)) %>%
      .[ , allele := alleles[[i]]]
  }) %>%
    bind_rows() %>%
    as_tibble()
  
  ### Generate NetMHCpan commands
  cmds <- peptide_exports_files %>%
    mutate(N_mer = params$N_mers) %>%
    mutate(output = str_replace_all(input, "netMHCpan_input", "netMHCpan_output")) %>%
    mutate(output = str_replace_all(output, ".tsv", ".txt")) %>%
    mutate(output = str_replace_all(output, "_ch_", paste0("_ch_", allele))) %>%
    mutate(cmds = paste(netMHCpan,
                        "-BA", "-inptype 1",
                        "-a", allele,
                        "-l", N_mer, 
                        "-p -f", input,
                        ">", output,
                        "-v")) %>%
    as_tibble()
  
  ### Run NetMHCpan
  if (!dir.exists(paste0(dir_DB_PTM_mz, "/netMHCpan_output"))) {
    dir.create(paste0(dir_DB_PTM_mz, "/netMHCpan_output"))
  }
  bettermc::mclapply(cmds$cmds, mc.cores = Ncpu, mc.cleanup=T, mc.preschedule=T, function(x){
    system(x, intern = T)
  }) %>%
    suppressMessages()
  
  ### Load NetMHCpan results
  # Parameters
  Affinity_threshold <- IC50_aggregation_table %>%
    dplyr::select(`MHC-I_alleles`, Affinity_threshold) %>%
    mutate(Affinity_threshold = as.numeric(Affinity_threshold)) %>%
    left_join(cmds, by = c("MHC-I_alleles" = "allele")) %>%
    unique()
  
  list_IC50 <- lapply(X = cmds$output, FUN = function(x){
    
    if (file.exists(x)) {
      # Read and filter netMHCpan output
      binders_df <- fread(x, nThread = Ncpu, sep = "^",
                          blank.lines.skip = TRUE, col.names = "V1", 
                          nrows = 52) 
      
      if (nrow(binders_df) > 42) {
        binders_df <- fread(x, nThread = Ncpu, sep = "^",
                            blank.lines.skip = TRUE, col.names = "V1",
                            skip = 52) %>%
          pull(V1) %>%
          str_replace_all(pattern = "[[:space:]]+", " ") %>%
          str_replace_all(pattern = "\\*", "") %>%
          str_split_fixed(pattern = " ", n = Inf) %>%
          as.data.table() %>%
          lazy_dt() %>%
          rename(MHC = V2, peptide = V3, `Aff(nM)` = V16) %>%
          select(MHC, peptide, `Aff(nM)`) %>%
          mutate(`Aff(nM)` = as.numeric(`Aff(nM)`),
                 length = nchar(peptide)) %>%
          filter((!is.na(`Aff(nM)`))) %>%
          mutate(peptide = str_remove_all(peptide, "X")) %>%
          as.data.table()
      } else {
        binders_df = data.table(MHC = character(),
                                peptide = character(),
                                `Aff(nM)` = numeric(),
                                length = numeric()
        )
        missing_allele <- Affinity_threshold %>%
          dplyr::filter(output == x) %>%
          pull(`MHC-I_alleles`) %>%
          print()
        print(paste0("Warning: empty input! ","Make sure that ", missing_allele, " is present in netMHCpan-4.1/data/allelenames"))
      }
      return(binders_df)
    } 
  }) %>%
    rbindlist() 
  cat(as.character(Sys.time()), " - ", "Predicting MHC-I affinity: Done", "\n")
  
  ### Add the columns for cutoffs
  cat(as.character(Sys.time()), " - ", "Format IC50 prediction: Start", "\n")
  list_IC50 <- list_IC50 %>%
    lazy_dt() %>%
    select(c("peptide", "MHC", "Aff(nM)")) %>%
    unique() %>%
    tidyr::pivot_wider(id_cols = c("peptide"),
                       names_from = c("MHC"),
                       values_from = c("Aff(nM)")) %>%
    
    as.data.table()
  cat(as.character(Sys.time()), " - ", "Format IC50 prediction: Done", "\n")
  cat(as.character(Sys.time()), " - ", "Add MHC-I affinities to peptide table: Start", "\n")
  
  # Remove previously used columns 
  for (a in seq_along(Affinity_cols)) {
    pep[[Affinity_cols[[a]]]] <- NULL
  }
  
  pep <- setkey(pep, peptide)
  pep <- merge(pep, list_IC50, all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE, by = "peptide") %>%
    unique()
  
  cat(as.character(Sys.time()), " - ", "Add MHC-I affinities to peptide table: Done", "\n")
  cat(as.character(Sys.time()), " - ", "Estimate IC50 cutoffs: Start", "\n")
  for (fj in seq_along(IC50_aggregation_table$Filename)) {
    colname_file <- paste0("Predicted_binder:", IC50_aggregation_table$Filename[[fj]], ":", 
                           IC50_aggregation_table$`MHC-I_alleles`[[fj]])
    allele <- paste0(IC50_aggregation_table$`MHC-I_alleles`[[fj]])
    IC50_threshold <- as.numeric(IC50_aggregation_table$Affinity_threshold[[fj]])
    cat(paste(colname_file, "|",IC50_threshold, "nM", "\n"))
    
    exists_MW.RT <- paste0("MW.RT.exists:", IC50_aggregation_table$Filename[[fj]])
    exists_MW.RT
    
    pep[, get("colname_file") := pep[[get("exists_MW.RT")]] & fifelse(.SD < IC50_threshold, T, F), 
        .SDcols = get("allele")]
  }
  setnames(pep, old = alleles, new = paste0("Aff(nM):", alleles))
  cat(as.character(Sys.time()), " - ", "Estimate IC50 cutoffs: Done", "\n")
} else {
  cat(as.character(Sys.time()), " - ", "No MHC-I affinity will be predicted", "\n")
} # End IC50 prediction

### Clean-up
suppressWarnings(rm(Affinity_threshold, cmds, IC50_aggregation_table, IC50_aggregation_table_j, list_IC50, 
                    mods, fixed_mods, mods_fixed_j, mods_fm, peptide_exports, tmp_1D, tmp_2D, RCs, 
                    MW.exists.j, MW.RT.exists.j, MS_mass_list, n_chunks, n_netMHCpan_chunks, predict_RT, 
                    PTM_list, PTM_mode, tolerance, RT_tolerance, a, fj, colname_file, IC50_threshold, allele))
