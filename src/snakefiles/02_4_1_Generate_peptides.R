### ---------------------------------------------- Generate PCP/PSP ----------------------------------------------
# description:  Generate unique PCP/PSP sequences
#               
# input:        1. Confirmed existence of Proteome_chunk.fasta
#               2. Parameters: index_length, peptide length, min_intervening_sequence length, output directory
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences .parquet
#               - Protein-peptide mapping .parquet
#               
# author:       YH, JL

### ---------------------------- (1) Proteome pre-processing --------------------------------------
cat(as.character(Sys.time()), " - ", "Proteome chunk extraction: Start \n")
Uniprot_regex <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"

# Extract Uniprot headers like in mmseqs2
if (FALSE %in% (names(dat) %in% proteome_index$desc)) {
  names(dat) <- ifelse(stri_detect_regex(names(dat), pattern = Uniprot_regex),
                       stri_extract(names(dat), regex = Uniprot_regex),
                       names(dat))
}

# I/L redundancy
if (replace_I_with_L == TRUE) {
  cat(as.character(Sys.time()), " - ", "Replace I with L: TRUE \n")
  dat <- chartr("I", "L", dat)
} else {
  cat(as.character(Sys.time()), " - ", "Replace I with L: FALSE \n")
}

# replace nul strings with X
nul <- as(as(raw(1), "XRaw"), "AAString")
dat <- chartr(nul, "X", dat)

# Split by terminal stop codons (if starts or ends with *)
dat[subseq(dat, 1, 1) == "*"] <- subseq(dat[subseq(dat, 1, 1) == "*"], start = 2)
dat[subseq(dat, width(dat), width(dat)) == "*"] <- subseq(dat[subseq(dat, width(dat), width(dat)) == "*"],
                                                          start = 1, end = width(dat[subseq(dat, width(dat), width(dat)) == "*"]) - 1)

# Filter by minimal length and re-order by similarity from clustering
dat <- dat[width(dat) >= min_protein_length]

# Separate long and short entries
short_entries <- dat[width(dat) <= max_protein_length]
long_entries <- dat[width(dat) > max_protein_length]

if (length(long_entries) > 0) {
  cat(as.character(Sys.time()), " - ", "Found sequences longer than max_length: splitting with overlap \n")
  
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
    right_join(select(proteome_index, desc)) %>%
    filter(!(is.na(value) | value == ""))
}
dat_split <- dat_split[names(dat_split) %in% long_entry_names$value] 
cat(as.character(Sys.time()), " - ", "Proteome chunk extraction: Done \n")
rm(proteome_index)
cat(as.character(Sys.time()), " - ", "Removed FASTA index", "\n")

### ---------------------------- (2) Cleavage peptides: PCP --------------------------------------
print(Sys.time())
dat_sort <- dat_split[width(dat_split) >= min(Nmers)]

if (enzyme_type == "PCP" | (grepl("cis-PSP", enzyme_type) == TRUE)) {
  
  # Pre-filter and sort input sequences
  cat(as.character(Sys.time()), " - ", "Computing PCP length: ", Nmers, " peptide generation\n")
  dat_sort <- dat_sort[order(width(dat_sort), decreasing = T),]
  
  split_chunks <- rep(1:(5*Ncpu), length.out=length(dat_sort))
  PCP <- split(dat_sort, split_chunks) %>%
    bettermc::mclapply(mc.preschedule = FALSE, 
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
  cat(as.character(Sys.time()), " - ", "Computing PCP: post-processing\n")
  PCP <- rbindlist(lapply(PCP, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  PCP <- PCP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, pattern = paste0("[^", str_c(AA, collapse = ""), "]")), .(peptide, protein), keyby=index] %>%
    .[, enzyme := "PCP"] %>%
    .[, length := str_length(peptide), keyby=index] 
  
  peptides <- PCP
  cat(as.character(Sys.time()), " - ", "Computed PCP\n")
} 

### ---------------------------- (3) Spliced peptides: PSP --------------------------------------
if (grepl("cis-PSP", enzyme_type) == TRUE) {
  cat(as.character(Sys.time()), " - ", "Computing cis-PSP: longest sequence index\n")
  
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
  cat(as.character(Sys.time()), " - ", "Computing cis-PSP: shorter sequences\n")
  
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
  cat(as.character(Sys.time()), " - ", "Computing cis-PSP: peptide generation\n")
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
  cat(as.character(Sys.time()), " - ", "Computing cis-PSP: post-processing\n")
  PSP <- setDT(stack(PSP))
  setnames(PSP, c("peptide", "protein"))
  PSP <- PSP %>% 
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, pattern = paste0("[^", str_c(AA, collapse = ""), "]")), .(peptide, protein), keyby=index] %>%
    .[, enzyme := "PSP", keyby=index] %>%
    .[, length := nchar(peptide), keyby=index]
  peptides <- rbindlist(list(PCP, PSP))
  
  cat(as.character(Sys.time()), " - ", "Computed PCP/PSP\n")
} 

### ---------------------------- (4) Cleave enzymatic digestions --------------------------------------
if (str_detect(enzyme_type, str_c(enzymes, collapse = "|"))) {
  
  enzym <- str_split_fixed(enzyme_type, "_", 2)[,1]
  missedCleavages <- str_split_fixed(enzyme_type, "_", 2)[,2]
  
  cat("Computing digestions by", enzym, "with", missedCleavages, "missed cleavages", "\n")
  cat("Peptides of length:", min(Nmers), "-", max(Nmers), "will be kept", "\n")
  
  # Generate in silico enzymatic digestions
  dat_sort <- dat[width(dat) >= min(Nmers)]
  
  ### Add the same sequences without N-terminal Met
  add_X <- dat_sort[as.character(subseq(dat_sort, 1, 1)) %in% c("M", "X")]
  add_X <- subseq(add_X, start = 2, end = ifelse(width(add_X) >= 2*max(Nmers), 2*max(Nmers), width(add_X)))
  dat_sort <- c(dat_sort, add_X)
  
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
  cat(as.character(Sys.time()), " - ", "Computing", enzyme_type, ": post processing\n")
  enzym_out <- rbindlist(lapply(enzym_out, as.data.table)) %>%
    lazy_dt() %>%
    select(value, group_name) %>%
    rename(protein = group_name,
           peptide = value) %>%
    as.data.table()
  
  enzym_out <- enzym_out %>%
    .[,index := str_sub(peptide, start = 1, end = index_length)] %>%
    .[!str_detect(peptide, pattern = paste0("[^", str_c(AA, collapse = ""), "]")), .(peptide, protein, index), keyby=index] %>%
    .[,length := str_length(peptide), keyby=index] %>%
    .[inrange(length, lower = min(Nmers), upper = max(Nmers), incbounds=TRUE), .(peptide, length, protein), keyby=index] %>%
    .[,enzyme := enzyme_type]
  
  peptides <- enzym_out
  if (nrow(enzym_out) == 0) {
    rm(enzym_out)
  }
  cat(as.character(Sys.time()), " - ", "Computed", enzyme_type, " \n")
}

### Format the protein-peptide table
cat(as.character(Sys.time()), " - ", "Formating the protein-peptide table", " \n")
peptides <- peptides[, protein := as.character(protein), keyby=index]
setkey(peptides, index, length, enzyme, peptide, protein)
peptides <- unique(peptides)

### Clean up
suppressWarnings(rm(dat, dat_sort, dat_split, add_X, index_list_result, inputs, PCP, PSP, enzym_out))

### ------------------------------------------ (5) Save peptide mapping ------------------------------------------
if (exists("peptides")) {
  cat(as.character(Sys.time()), " - ", "Saving peptide mapping \n")
  
  peptides <- peptides[, c("proteome", "MiSl", "chunk") := list(chunk_params$proteome, 
                                                                chunk_params$MiSl, 
                                                                chunk_params$chunk)] 
  peptides[1:nrow(peptides), 
           write_dataset(group_by(as_arrow_table(.SD), 
                                  length, proteome, enzyme, MiSl, chunk), 
                         path = paste0(dir_DB_exhaustive,  "/peptide_mapping/index=", .BY), 
                         existing_data_behavior = "overwrite",
                         format = "parquet", 
                         max_partitions = 10240L,
                         use_dictionary = FALSE,
                         compression = "lz4"), 
           by=index, 
           .SDcols=colnames(peptides)]
  
  cat(as.character(Sys.time()), " - ", "Saved peptide mapping \n")
}

### Clean up
gc()
