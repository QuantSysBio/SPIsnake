### ---------------------------------------------- SPIsnake: main ----------------------------------------------
# description:  Create unique .FASTA outputs
#               
# input:        1. Generate_peptides rules are done
#               2. Peptide sequences are saved as .FASTA and .CSV per chunk. 
#                  For PTM-modified peptides, there should exist arrow datasets
#               3. Experiment_design, Master_table_expanded
# output:       
#               - Unique peptides per Biological group as defined in Experiment_design
#               
# author:       YH, JL, HPR

### ---------------------------- Common variables ----------------------------
monoisotopic_masses <- data.frame(
  AA = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"),
  mass = c(71.037114,103.009185,115.026943,129.042593,147.068414,57.021464,137.058912,113.084064,128.094963,113.084064,131.040485,114.042927,97.052764,128.058578,156.101111,87.032028,101.047679,99.068414,186.079313,163.06332)
)

# Proteinogenic AAs
AA = monoisotopic_masses$AA

# Cleaver params
enzymes <- c("arg-c proteinase", "asp-n endopeptidase", "bnps-skatole", "caspase1", "caspase2", "caspase3", "caspase4", "caspase5", "caspase6", "caspase7", "caspase8", "caspase9", "caspase10", "chymotrypsin-high", "chymotrypsin-low", "clostripain", "cnbr", "enterokinase", "factor xa", "formic acid", "glutamyl endopeptidase", "granzyme-b", "hydroxylamine", "iodosobenzoic acid", "lysc", "lysn", "neutrophil elastase", "ntcb", "pepsin1.3", "pepsin", "proline endopeptidase", "proteinase k", "staphylococcal peptidase i", "thermolysin", "thrombin", "trypsin")

### ---------------------------- Proteome pre-processing ----------------------------
Split_max_length <- function(x, max_length=500, overlap_length=MiSl*2, chunk_positions = F){
  lapply(x, function(x){
    
    if (nchar(x) <= max_length) {
      protein_chunks <- x
    } else {
      # Extract coordinates
      totalsize = 200*nchar(x) 
      
      at <- breakInChunks(totalsize = totalsize, chunksize = max_length)
      at <- successiveIRanges(width = width(at), gapwidth = -overlap_length)
      
      # Subset last relevant IRanges
      keep <- which((at@start + at@width) <= nchar(x))
      keep <- c(keep, max(keep) + 1)
      at <- at[keep]
      
      ### Correct the last coordinate
      at@width[length(at@width)] <- as.integer(nchar(x) - at@start[length(at@start)] + 1)
      
      # Extract sequence
      protein_chunks <- extractAt(x, at)
      
      # Update names
      if (chunk_positions) {
        names(protein_chunks) <- paste0("|chunk:", at@start, "-", at@start + at@width - 1)
      }
    }
    return(protein_chunks)
  })
}

Save_prot_chunk_biostrings <- function(dat, 
                                       nF=nF, 
                                       Pi=Pi,
                                       MiSl=MiSl,
                                       orderedProteomeEntries=orderedProteomeEntries, 
                                       directory=directory,
                                       proteome_name=proteome_name,
                                       protein_counter_start=protein_counter_start,
                                       protein_counter_end=protein_counter_end,
                                       maxE=maxE){
  for(j in 1:nF){
    start=Pi[j]+1
    end=Pi[j+1]
    
    writeXStringSet(x = dat[orderedProteomeEntries[start:end]], 
                    filepath = paste0(directory, "/", "maxE_", maxE, "_MiSl_", MiSl, "_", proteome_name, "_", protein_counter_start, "_", protein_counter_end, "_", j,".fasta"))
  }
}

### ---------------------------- PCP/PSP generation ----------------------------
seq_vectorized <- Vectorize(seq.default, vectorize.args = c("from", "to"), SIMPLIFY = F) 

# compute all PCP with length <=nmer
computeCPomplete <- function(L,nmer){
  maxL = nmer+1
  CP = numeric()
  
  for(i in 1:L){
    CP = rbind(CP, cbind(rep(i,length(c(i:min(L,(i+maxL-2))))),
                         c(i:min(L,(i+maxL-2)))))
  }
  return(CP)
  rm(maxL, CP)
}

# translate PCP
translateCP <- function(CP,peptide){
  
  CPseq = rep(NA,dim(CP)[1])
  for(i in 1:dim(CP)[1]){
    CPseq[i] = paste(peptide[CP[i,1]:CP[i,2]],sep="",collapse="")
  }
  return(CPseq)
  rm(CPseq)
}

# compute all PSP with length == nmer
computeSPcomplete <- function(cp,maxL,minL,MiSl){
  SP = numeric()
  N = dim(cp)[1]
  NN = 5 * (10**6)
  
  SP = matrix(NA,NN,4)
  a = 1
  # repeat as many times as you have cp
  for(i in 1:N){
    temp1 = rep(cp[i,1],N)
    temp2 = rep(cp[i,2],N)
    temp3 = cp[,1]
    temp4 = cp[,2]
    
    L = temp4-temp3+temp2-temp1+2
    ind = which(((temp3-temp2)==1)|(L>maxL)|(L<minL)|((temp3-temp2)>(MiSl+1))|((temp1-temp4)>(MiSl+1))|((temp3<=temp2)&(temp4>=temp1)))
    
    if(length(ind)>0){
      temp1 = temp1[-ind]
      temp2 = temp2[-ind]
      temp3 = temp3[-ind]
      temp4 = temp4[-ind]
    }
    
    #
    if(length(temp1)!=1){
      if((a+length(temp1)-1)>NN){
        ## this looks slow but is rarely called so it is fine
        SP = rbind(SP,matrix(NA,NN,4))
      }
      
      if(length(temp1)>0){
        SP[c(a:(a+length(temp1)-1)),1] = temp1
        SP[c(a:(a+length(temp1)-1)),2] = temp2
        SP[c(a:(a+length(temp1)-1)),3] = temp3
        SP[c(a:(a+length(temp1)-1)),4] = temp4
      }
      a = a+length(temp1)
    }
  }
  # remove all empty lines from SP
  if(a<NN){
    SP = SP[-c(a:NN),]
  }
  return(SP)
  rm(temp, N, NN, L, ind, SP, a)
}

CutAndPaste_seq_return_sp <- function(inputSequence,nmer,MiSl){
  ### Makes PCP, cis-PSP, revcis-PCP from a given input
  
  results = list()
  peptide = strsplit(inputSequence,"")[[1]]
  L = length(peptide)
  
  if(L>nmer){
    
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L, nmer)
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1) ==nmer)
    cpNmer = cp[index,]
    
    if(length(cp[,1])>1){
      
      ## these are the indicese
      sp = computeSPcomplete(cp,maxL=nmer,minL=nmer,MiSl=MiSl)
      return(sp)
    }
  }
}

CutAndPaste_seq_return_sp_vec <- Vectorize(CutAndPaste_seq_return_sp, 
                                           vectorize.args = "inputSequence", 
                                           SIMPLIFY = F, 
                                           USE.NAMES = F)

Generate_PSP_2 <- function(protein_inputs, Nmers){
  stri_join(extractAt(protein_inputs[[1]], IRanges(start = protein_inputs[[2]][,1], 
                                                   end =   protein_inputs[[2]][,2])), 
            extractAt(protein_inputs[[1]], IRanges(start = protein_inputs[[2]][,3], 
                                                   end =   protein_inputs[[2]][,4])))
}

### ---------------------------- MW and PTMs ----------------------------

read_MW_file <- function(file, num_threads){
  fread(file = file, 
        sep = " ", nThread = num_threads, data.table = T,
        col.names = c("Precursor_mass","RT")) %>%
    lazy_dt() %>%
    mutate(MW_Min = Precursor_mass - Precursor_mass * tolerance * 10 ** (-6)) %>%
    mutate(MW_Max = Precursor_mass + Precursor_mass * tolerance * 10 ** (-6)) %>%
    ### Min/sec for RT
    # mutate(RT = RT / 60) %>%
    mutate(RT_Min = RT - RT_tolerance) %>%
    mutate(RT_Max = RT + RT_tolerance) %>%
    select(-Precursor_mass)  %>%
    unique() %>%
    as.data.table()
}

computeMZ_biostrings <- function(seq, AAs = monoisotopic_masses$AA, masses = monoisotopic_masses$mass){
  seq = AAStringSet(seq)
  
  if(length(seq) < 1){
    MW = NA
  }
  if(length(seq) > 0){
    aa3 = letterFrequency(seq, letters = AAs) %*% diag(masses)
    MW=rowSums(aa3)+18.01056
    ## add NA for zero values
    # MW=round(MW, digits = 5)
    MW[MW==18.01056] = NA
  }
  return(MW)
}

### ---------------------------- RT filtering ----------------------------
regression_stats <- function(obs, pred, quantile_user = 0.99){
  # lm
  data = data.frame(pred = pred, 
                    obs = obs)
  pred.lm = lm(pred ~ obs, data = data)
  
  # metrics
  # correlation coefficients
  pc = cor(obs, pred, method = "pearson")
  sm = cor(obs, pred, method = "spearman")
  
  # mean squared error
  mse =  round(mean((obs - pred)^2), 4)
  # root mean squared error
  rmse = round(sqrt(mse), 4)
  # mean absolute deviation
  mae =  round(mean(abs((obs - pred))), 4)
  # q95
  q95 =  round(quantile(abs((obs - pred)), probs = 0.95), 4)
  # q99
  q99 =  round(quantile(abs((obs - pred)), probs = 0.99), 4)
  # quantile cutoff by user
  quantile_user =  round(quantile(abs((obs - pred)), probs = quantile_user), 4)
  
  # sumarize
  all.metrics = c(summary(pred.lm)$r.squared, pc, mse, rmse, mae, q95, q99, quantile_user)
  names(all.metrics) = c("Rsquared", "PCC", "MSE", "RMSE", "MAE", "quantile_95", "quantile_99", "quantile_user")
  
  return(all.metrics)
}


### ---------------------------- PTM generation ----------------------------
find_ka <- Vectorize(function(aa, mods_input){
  # which(mods_input$Site %in% aa & mods_input$Position=="Anywhere")
  which(mods_input$Site %chin% aa & mods_input$Position=="Anywhere")
}, SIMPLIFY = T, USE.NAMES = F)

getPTMcombinations_fast <- function(s = peptide, m = MW, NmaxMod=max_variable_PTM, mods_input=mods){
  
  modIndex = c(which(mods_input$Site=="N-term" & mods_input$Position=="Any N-term"),
               which(mods_input$Site=="C-term" & mods_input$Position=="Any C-term"),
               unlist(find_ka(aa = strsplit(s, split="")[[1]], 
                              mods_input = list(mods_input))))
  
  modId <- mods_input$Id[modIndex]
  modDelta <- as.numeric(mods_input$MonoMass[modIndex])
  
  # get all combinations of NmaxMod modifications
  combi = list()
  combi[[1]] = c(1:length(modIndex))
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod, length(modIndex))){
      combi[[i]] <- arrangements::combinations(c(1:length(modIndex)), i)
    }
  }
  
  deltaMass <- list()
  IDs <- list()
  deltaMass[[1]] <- modDelta[combi[[1]]]
  IDs[[1]] <- modId[combi[[1]]]
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod,length(modIndex))){
      deltaMass[[i]] <- rowSums(matrix(data = modDelta[combi[[i]]],
                                       nrow = dim(combi[[i]])[1], 
                                       ncol = i))
      
      IDs[[i]] <- stri_join_list(x = split(modId[combi[[i]]], f = rep(1:dim(combi[[i]])[1], times = i)), sep=";")
    }
  }
  
  # generate final mod sequences with delta Masses
  PTMcombinations <- data.table(peptide = s,
                                ids = unlist(IDs),
                                MW = m + unlist(deltaMass))
  return(PTMcombinations)
}
getPTMcombinations_fast_vec <- Vectorize(getPTMcombinations_fast, vectorize.args = c("s", "m"), SIMPLIFY = F)

getPTMcombinations_fixed <- function(s = peptide, m = MW, mods_input=mods){
  
  modIndex <- c(which(mods_input$Site=="N-term" & mods_input$Position=="Any N-term"),
                which(mods_input$Site=="C-term" & mods_input$Position=="Any C-term"),
                unlist(find_ka(aa = strsplit(s, split="")[[1]], 
                               mods_input = list(mods_input))))
  
  # generate final mod sequences with delta Masses
  PTM_fixed <- data.table(peptide = s,
                          ids = str_c(sort(na.omit(mods_input$Id[modIndex])), collapse = ";"),
                          MW = m + sum(na.omit(as.numeric(mods_input$MonoMass[modIndex]))))
  return(PTM_fixed)
}
getPTMcombinations_fixed_vec <- Vectorize(getPTMcombinations_fixed, vectorize.args = c("s", "m"), SIMPLIFY = F)

getPTMcombinations_fixed_mass <- function(s = peptide, m = MW, mods_input=mods){
  
  modIndex <- c(which(mods_input$Site=="N-term" & mods_input$Position=="Any N-term"),
                which(mods_input$Site=="C-term" & mods_input$Position=="Any C-term"),
                unlist(find_ka(aa = strsplit(s, split="")[[1]], 
                               mods_input = list(mods_input))))
  
  # final peptide mass
  MW <- m + sum(na.omit(as.numeric(mods_input$MonoMass[modIndex])))
  return(MW)
}
getPTMcombinations_fixed_mass_vec <- Vectorize(getPTMcombinations_fixed_mass, vectorize.args = c("s", "m"), SIMPLIFY = F)


# -------------------------------------- Arrow aggregation --------------------------------------
try_with_timeout <- function(expr, timeout = timeout, ...){
  require(R.utils)
  tryCatch({withTimeout({expr}, timeout = timeout)
  }, TimeoutException = function(ex) {
    message(paste0("Timeout limit: ", 10*timeout, " seconds"))
  },
  error=function(cond) {
    message("Error")
    return(NULL)
  })
}

try_with_timeout_restart <- function(expr, retries, timeout, wait_on_restart=0) {
  require(R.utils)
  withRestarts(expr = tryCatch({
    res <- withTimeout({expr}, timeout = timeout)
  }, TimeoutException = function(ex) {
    message(paste0("Timeout limit: ", timeout*10, " seconds"))
    Sys.sleep(wait_on_restart)
    invokeRestart("retry")
  },
  error=function(e) { 
    Sys.sleep(wait_on_restart)
    invokeRestart("retry")
  }),
  retry = function() { 
    message(paste0("Retrying. Attepmts left : "), retries)
    stopifnot(retries > 0)
    try_with_timeout_restart(expr = expr, 
                             retries = retries-1, 
                             timeout = timeout, 
                             wait_on_restart = wait_on_restart)
  })
}
memfree <- function(){
  as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) / 2^20
}

terminate_duckdb <- function(conn = conn){
  cat(as.character(Sys.time()), " - ", "disconnecting and shutting down duckdb ", "\n")
  duckdb::duckdb_shutdown(duckdb::duckdb())
  rm(conn, DB_duck)
}

clean_duckdb_RAM <- function(){
  RAM_free <- memfree()
  if (RAM_free / Max_RAM < 0.2) {
    cat(as.character(Sys.time()), " - ", "Starting clean up:", round(RAM_free),"Gb left free in RAM\n")
    if (exists("conn")) {
      terminate_duckdb()
    }
    gc(full = T, reset = T)
    cat(as.character(Sys.time()), " - ", "Done clean up", "\n")
    
    conn <- DBI::dbConnect(duckdb(), dbname = paste0("arrow_unique.duckdb.", aggregation_batch_i),
                           config=list("memory_limit"= paste0(duckdb_RAM, "GB"),
                                       "temp_directory" = duckdb_temp_dir))
    cat(as.character(Sys.time()), " - ", "Connected to the database", "\n")
  }
}

get_distict_duckdb <- function(input_path = paste0(dir_DB_PTM_mz, "peptide_seqences/"),
                               DB_groups_i = DB_groups_i,
                               table_name = table_name,
                               conn = conn
){
  paste0(input_path, "/index=", DB_groups_i$index, "/length=", DB_groups_i$length) %>%
    open_dataset() %>%
    select(-chunk) %>%    
    to_duckdb(con = conn, 
              table_name = table_name) %>%
    mutate(index = local(DB_groups_i$index[[1]])) %>%
    mutate(length = local(DB_groups_i$length[[1]])) %>%
    window_order() %>%
    group_by(proteome, enzyme, MiSl) %>%
    distinct(peptide, .keep_all = T) %>%
    mutate(index = local(DB_groups_i$index[[1]])) %>%
    window_order() %>%
    group_by(index, length, proteome, enzyme, MiSl) %>%
    collapse() 
}

try_with_time_limit <- function(expr, cpu = Inf, elapsed = Inf) {
  y <- try({setTimeLimit(cpu, elapsed); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) NULL else y 
}

aggregate_peptide_mapping <- function(input_path = paste0(dir_DB_exhaustive, "/peptide_mapping/"), 
                                      DB_groups_i,
                                      output_dir = "peptide_mapping"){
  paste0(input_path, "/index=", DB_groups_i$index, "/length=", DB_groups_i$length) %>%
    open_dataset() %>%
    select(-chunk) %>%  
    to_duckdb() %>%
    mutate(index = local(DB_groups_i$index[[1]])) %>%
    mutate(length = local(DB_groups_i$length[[1]])) %>%
    to_arrow() %>%
    group_by(index, length, proteome, enzyme, MiSl) %>%
    write_dataset(path = output_dir,
                  existing_data_behavior = "overwrite",
                  format = "parquet",
                  max_partitions = 10240L,
                  max_rows_per_file = as.integer(2 * 10^8),
                  compression = "lz4") 
}

make_fasta <- function(DB_arrow_files = DB_arrow_files,
                       groups_j = groups_j,
                       j = j, 
                       fasta_Biological_group = "",
                       filter_prefix = filter_prefix,
                       filter_suffix = filter_suffix,
                       path_out = ""){
  
  out <- DB_arrow_files %>%
    right_join(groups_j[j,], by = c("index", "length")) %>%
    mutate(dataset_path = str_split_fixed(filename, pattern = "/proteome=", n = 2)[,1]) %>%
    pull(dataset_path) %>%
    unique() %>%
    open_dataset()
  cat(as.character(Sys.time()), ":", "Opened dataset", "\n")
  
  # Data-driven filters (optional)
  if (!is.null(filter_prefix) & !is.null(filter_suffix)) {
    out <- out %>%
      mutate(index = groups_j[j,]$index, 
             length = groups_j[j,]$length) %>%
      group_by(index, length, proteome, enzyme, MiSl) %>%
      dplyr::select(peptide, index, length, proteome, enzyme, MiSl,
                    starts_with(paste0(filter_prefix, filter_suffix))) %>%
      filter(if_any(starts_with(paste0(filter_prefix, filter_suffix)), ~ .)) %>%
      collapse() %>%
      group_by(proteome, enzyme, MiSl) %>%
      select(peptide, proteome, enzyme, MiSl, index, length) %>%
      collapse() %>%
      compute()
    cat(as.character(Sys.time()), ":", "Done filtering", "\n")
  } else {
    out <- out %>%
      mutate(index = groups_j[j,]$index, 
             length = groups_j[j,]$length) %>%
      group_by(proteome, enzyme, MiSl) %>%
      select(peptide, proteome, enzyme, MiSl, index, length) %>%
      collapse() %>%
      compute()
    cat(as.character(Sys.time()), ":", "No filtering", "\n")
  }
  
  # Save FASTA
  out <- out %>%
    to_duckdb() %>%
    mutate(peptide = paste0(">", index, "_", length, "_", row_number(), "\n", peptide)) %>%
    mutate(group = paste0(enzyme, "_", MiSl, "_", proteome)) %>%
    ungroup() %>%
    select(peptide, group) %>%
    collapse() %>%
    compute() 
  cat(as.character(Sys.time()), ":", "Formatted FASTA", "\n")
  
  out <- out %>%
    collect() %>%
    setDT() 
  print(head(out))
  cat(as.character(Sys.time()), ":", "Loaded data.table", "\n")
  
  if (nrow(out) > 0) {
    fasta_Biological_group <- ifelse(nchar(fasta_Biological_group) > 0, 
                                     yes = paste0("_", fasta_Biological_group),
                                     no = "")
    cat(as.character(Sys.time()), ":", "Writing FASTA", "\n")
    
    out[, fwrite(.SD, paste0(path_out, "_", group, fasta_Biological_group, ".fasta"),
                 sep = "\n", col.names = FALSE, append = T, 
                 quote = FALSE, nThread = Ncpu),
        by=list(group), .SDcols=c("peptide")]
  }
  return(TRUE)
}

make_fasta_batched <- function(DB_arrow_files = DB_arrow_files,
                               groups = groups,
                               groups_j = groups_j,
                               fasta_Biological_group = "",
                               filter_prefix = NULL,
                               filter_suffix = NULL,
                               timeout = 60,
                               retries = 20,
                               wait_on_restart = wait_on_restart,
                               path_out = ""){
  
  # Remove FASTA files with matching names
  suppressWarnings(file.remove(paste0(path_out, "_",
                                      groups$enzyme, "_",
                                      groups$MiSl, "_", 
                                      groups$proteome, 
                                      ifelse(nchar(fasta_Biological_group) > 0, paste0("_", fasta_Biological_group),""),
                                      ".fasta")))
  
  for (j in seq_along(groups_j$length)) {
    
    cat(as.character(Sys.time()), ":", groups_j$index[[j]], groups_j$length[[j]], "|", j, "/", length(groups_j$length), "\n")
    out <- try_with_timeout_restart({
      make_fasta(DB_arrow_files = DB_arrow_files, 
                 j = j, 
                 groups_j = groups_j,
                 path_out = path_out,
                 fasta_Biological_group = fasta_Biological_group,
                 filter_prefix = filter_prefix,
                 filter_suffix = filter_suffix)
    }, 
    wait_on_restart = wait_on_restart,
    timeout = timeout, 
    retries = retries)
  } # End index-length
  return(out)
}

### ---------------------------- Logging ----------------------------
SPIsnake_log <- function(memory_profile = F){
  print("----- memory usage by Slurm -----")
  jobid = system("echo $SLURM_JOB_ID")
  system(paste0("sstat ", jobid)) %>%
    print()
  
  system("sacct --format='JobID,JobName,State,Elapsed,AllocNodes,NCPUS,NodeList,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask,ReqMem,MaxDiskWrite'") %>%
    print()
  
  if (memory_profile == T) {
    print("----- memory usage by R -----")
    memory.profile() %>%
      print()
  }
  
  print("----- connections -----")
  showConnections() %>%
    print()
  
  print("----- removing cluster -----")
  if (exists("cl")) {
    print(cl)
    parallel::stopCluster(cl)
  }
  
  print(sessionInfo())
  if (exists("snakemake")) {
    sink()
  }
} 

