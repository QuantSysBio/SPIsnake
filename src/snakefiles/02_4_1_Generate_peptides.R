### ---------------------------------------------- Generate PCP/PSP ----------------------------------------------
# description:  Generate unique PCP/PSP sequences
#               
# input:        1. Confirmed existance of Proteome_chunk.fasta
#               2. Parameters: index_length, peptide length, min_intervening_sequence length, output directory
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences .fst
#               - Protein-peptide mapping .fst
#               - Peptide generation stats .csv
#               
# author:       YH, JL, KP

### ---------------------------- (1) Proteome pre-processing --------------------------------------
cat("Proteome chunk extraction: Start \n", as.character(Sys.time()), " \n")
# Extract Uniprot headers like in mmseqs2
if (FALSE %in% (names(dat) %in% proteome_index$desc)) {
  names(dat) <- ifelse(stri_detect_regex(names(dat), pattern = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                       stri_extract(names(dat), regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"),
                       names(dat))
}

# I/L redundancy
if (replace_I_with_L == TRUE) {
  cat("Replace I with L: TRUE \n", as.character(Sys.time()), "\n")
  dat <- chartr("I", "L", dat)
} else {
  cat("Replace I with L: FALSE \n", as.character(Sys.time()), "\n")
}

# Filter by minimal length and re-order by similarity from clustering
dat <- dat[width(dat) >= min_protein_length]

# Separate long and short entries
short_entries <- dat[width(dat) <= max_protein_length]
long_entries <- dat[width(dat) > max_protein_length]

if (length(long_entries) > 0) {
  cat("Found sequences longer than max_length: splitting with overlap \n", as.character(Sys.time()), "\n")
  
  # Split long entries into chunks
  long_entries <- split(long_entries, f = rep(1:Ncpu, length.out=length(long_entries), each=ceiling(length(long_entries)/Ncpu)))
  if (Ncpu > 1) {
    long_entries <- bettermc::mclapply(X = long_entries,
                                       FUN = Split_max_length,
                                       max_length = max_protein_length,
                                       overlap_length = 2*params$Max_Interv_length,
                                       mc.preschedule = T,
                                       mc.cleanup = T,
                                       mc.retry = 3,
                                       mc.cores = Ncpu)
  } else {
    long_entries <- lapply(long_entries, Split_max_length, max_length=max_protein_length, overlap_length=2*params$Max_Interv_length)
  }
  long_entries <- lapply(long_entries, function(x){unlist(AAStringSetList(x), use.names = T)})
  long_entries <- unlist(AAStringSetList(long_entries), use.names = F)
  
  # Re-merge and order
  dat_split <- c(short_entries, long_entries)
  long_entry_names <- names(dat_split) %>%
    as_tibble() %>%
    mutate(desc = str_split_fixed(value, "\\.\\|chunk", 2)[,1]) %>%
    right_join(select(proteome_index, desc)) %>%
    filter(!(is.na(value) | value == ""))
} else {
  dat_split <- dat
  long_entry_names <- names(dat_split) %>%
    as_tibble() %>%
    mutate(desc = value) %>%
    right_join(select(proteome_index, desc))
}
dat_split <- dat_split[long_entry_names$value]
cat("Proteome chunk extraction: Done\n", as.character(Sys.time()), "\n")

### ---------------------------- (2) Cleavage peptides: PCP --------------------------------------
print(Sys.time())
dat_sort <- dat_split[width(dat_split) >= min(Nmers)]

if (enzyme_type == "PCP" | (grepl("cis-PSP", enzyme_type) == TRUE)) {
  
  # Pre-filter and sort input sequences
  cat("Computing PCP length: ", Nmers, " peptide generation\n", as.character(Sys.time()), "\n")
  dat_sort <- dat_sort[order(width(dat_sort), decreasing = T),]
  
  split_chunks <- rep(1:(5*Ncpu), length.out=length(dat_sort))
  PCP <- split(dat_sort, split_chunks) %>%
    bettermc::mclapply(mc.preschedule = TRUE, 
                       mc.cores = Ncpu, 
                       mc.cleanup = TRUE, 
                       mc.force.fork = TRUE, 
                       mc.retry = 3, 
                       Nmers = Nmers,
                       FUN = function(x, Nmers){
                         extractAt(x, IRangesList(start = seq_vectorized(from = 1, to = width(x)-Nmers+1),
                                                  end = seq_vectorized(from = Nmers, to = width(x))))
                       })
  # Tidy format
  cat("Computing PCP: post-processing\n", as.character(Sys.time()), " \n")
  PCP <- rbindlist(lapply(PCP, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  setkey(PCP, peptide)
  PCP <- PCP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein), by=index] %>%
    .[, enzyme := "PCP"] %>%
    .[, length := str_length(peptide), by=index] 

  peptides <- PCP
  cat("Computed PCP\n", as.character(Sys.time()), " \n")
} 

### ---------------------------- (3) Spliced peptides: PSP --------------------------------------
if (grepl("cis-PSP", enzyme_type) == TRUE) {
  cat("Computing cis-PSP: longest sequence index\n", as.character(Sys.time()), " \n")
  
  
  ### Generate indices
  longest_index <- parallel::mcmapply(FUN = CutAndPaste_seq_return_sp_vec, 
                                      inputSequence = paste0(sample(c("A","T","G","C"), max(width(dat_sort)), replace=T), collapse=""), 
                                      nmer = Nmers, 
                                      MiSl = params$Max_Interv_length, 
                                      SIMPLIFY = F, USE.NAMES = F,
                                      mc.cores = Ncpu, 
                                      mc.preschedule = T) %>%
    unlist(recursive = F) %>%
    as.data.table()
  
  ### Match indices to proteins of corresponding length
  cat("Computing cis-PSP: shorter sequences\n", as.character(Sys.time()), " \n")
  inputs <- lapply(seq_along(1:length(dat_sort)), function(x){
    out <- list()
    out[[1]] <- dat_sort[[x]]
    out[[2]] <- longest_index %>%
      lazy_dt() %>%
      filter(V1 <= width(dat_sort[x]) & V2 <= width(dat_sort[x]) & V3 <= width(dat_sort[x]) & V4 <= width(dat_sort[x])) %>%
      as.data.table() %>%
      as.matrix()
    return(out)
  })
  
  # Generate PSP sequences
  cat("Computing cis-PSP: peptide generation\n", as.character(Sys.time()), " \n")
  PSP <- bettermc::mclapply(X = inputs, 
                            Nmers = Nmers,
                            FUN = Generate_PSP_2,
                            mc.cores = Ncpu,
                            mc.cleanup = TRUE, 
                            mc.preschedule = TRUE, 
                            mc.force.fork = TRUE, 
                            mc.retry = 3)
  names(PSP) <- names(dat_sort)
  
  # Post-processing
  cat("Computing cis-PSP: post-processing\n", as.character(Sys.time()), " \n")
  PSP <- setDT(stack(PSP))
  setnames(PSP, c("peptide", "protein"))
  setkey(PSP, peptide)
  PSP <- PSP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein), by=index] %>%
    .[, enzyme := "PSP"] %>%
    .[, length := str_length(peptide), by=index]
  peptides <- rbindlist(list(PCP, PSP))
  
  cat("Computed PCP/PSP\n", as.character(Sys.time()), " \n")
} 

### ---------------------------- (4) Cleave enzymatic digestions --------------------------------------
if (str_detect(enzyme_type, str_c(enzymes, collapse = "|"))) {
  
  enzym <- str_split_fixed(enzyme_type, "_", 2)[,1]
  missedCleavages <- str_split_fixed(enzyme_type, "_", 2)[,2]
  
  
  cat("Computing digestions by", enzym, "with", missedCleavages, "missed cleavages")
  cat("Peptides of length:", min(Nmers), "-", max(Nmers), "will be kept")
  
  # Generate in silico enzymatic digestions
  dat_sort <- dat[width(dat) >= min(Nmers)]
  
  ### Add the same sequences without N-terminal Met
  {
    add_X <- dat_sort[str_starts(dat_sort, "M|X")]
    add_X <- subseq(add_X, start = 2, end = ifelse(width(add_X) >= 2*max(Nmers), 2*max(Nmers), width(add_X)))
    dat_sort <- c(dat_sort, add_X)
  }
  
  split_chunks <- rep(1:(5*Ncpu), length.out=length(dat_sort))
  if (enzym == "trypsin") {
    enzym_out <- split(dat_sort, split_chunks) %>%
      bettermc::mclapply(FUN = cleave,
                         missedCleavages = 0:missedCleavages,
                         custom = custom_trypsin, 
                         unique = FALSE,
                         mc.cores = Ncpu,
                         mc.cleanup = TRUE, 
                         mc.preschedule = TRUE, 
                         mc.force.fork = TRUE, 
                         mc.retry = 3)
  } else {
    enzym_out <- split(dat_sort, split_chunks) %>%
      bettermc::mclapply(FUN = cleave,
                         enzym = enzym, 
                         missedCleavages = 0:missedCleavages,
                         unique = FALSE,
                         mc.cores = Ncpu,
                         mc.cleanup = TRUE, 
                         mc.preschedule = TRUE, 
                         mc.force.fork = TRUE, 
                         mc.retry = 3)
  }
  
  # Tidy format
  cat("Computing", enzyme_type, ": post processing\n", as.character(Sys.time()), " \n")
  enzym_out <- rbindlist(lapply(enzym_out, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  setkey(enzym_out, peptide)
  enzym_out <- enzym_out %>% 
    unique() %>%
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, exclusion_pattern), .(peptide, protein, index), by=index] %>%
    .[,length := str_length(peptide), by=index] %>%
    .[inrange(length, lower = min(Nmers), upper = max(Nmers), incbounds=TRUE), .(peptide, length, protein), by=index] %>%
    .[, enzyme := enzyme_type]
  
  peptides <- enzym_out
  if (nrow(enzym_out) == 0) {
    rm(enzym_out)
  }
  cat("Computed", enzyme_type, " \n", as.character(Sys.time()), " \n")
}

### Format the protein-peptide table
cat("Formating the protein-peptide table", " \n", as.character(Sys.time()), " \n")
peptides <- peptides[, protein := as.character(protein), by=index]
setkey(peptides, enzyme, index, length, peptide, protein)

### Clean up
suppressWarnings(rm(dat, dat_sort, dat_split, add_X, index_list_result, inputs, PCP, PSP, enzym_out))

### ------------------------------------------ (5) Save peptide mapping ------------------------------------------
if (exists("peptides")) {
  cat("Saving peptide mapping \n", as.character(Sys.time()), " \n")
  
  peptides[1:nrow(peptides)] %>%
    as_arrow_table() %>%
    mutate(proteome = chunk_params$proteome) %>%
    mutate(MiSl = chunk_params$MiSl) %>%
    mutate(chunk = chunk_params$chunk) %>%
    group_by(enzyme, MiSl, proteome, index, length, chunk) %>%
    write_dataset(path = paste0(dir_DB_exhaustive, "/peptide_mapping/"),
                  existing_data_behavior = "overwrite",
                  format = "parquet",
                  compression = "lz4")
  
  cat("Saved peptide mapping \n", as.character(Sys.time()), " \n")
}

### Clean up
gc()
