### ---------------------------------------------- Aggregate IC50  ----------------------------------------------
# description:  Join protein-peptide mapping, MW- RT- and binding affinity predictions.
#               Report search space sizes.
#               
# input:        1. metadata: Master_table_expanded, Experiment_design, cmd_netMHCpan
#               2. Candidate peptide binders filtered from netMHCpan output
#               3. Peptides with MW & RT information after 2D filter
#               
# output:       
#               - FASTA file per Biological group
#               - peptides and their properties per Biological group in .csv
#               - DB size at each filtering step (table and plot)
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log, split = TRUE)

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))
print(sessionInfo())

# {
#   ### setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
#   ### setwd("/data/SPIsnake/K562_expressed_PCP_SPIsnakeVersion20220218")
#   Master_table_expanded <- fread("results/DB_exhaustive/Master_table_expanded.csv")
#   Experiment_design <- fread("data/Experiment_design.csv")
#   dir_DB_exhaustive = "results/DB_exhaustive/"
#   dir_DB_PTM_mz = "results/DB_PTM_mz"
#   dir_IC50 = "results/IC50/"
#   dir_DB_out = "results/DB_out"
#   cmd_netMHCpan <- fread(paste0(dir_IC50, "cmd_netMHCpan.csv"))
#   fst_compression = 100
#   minimal_output_headers = TRUE
# }

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
gg <- list()

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

# fst compression level
fst_compression = as.integer(snakemake@params[["fst_compression"]])

# Header size
minimal_output_headers = as.logical(snakemake@params[["minimal_output_headers"]])

### ---------------------------- (1) Read inputs ----------------------------
# Experiment_design
Experiment_design <- fread(snakemake@input[["Experiment_design"]]) 
Master_table_expanded <- fread(snakemake@input[["Master_table_expanded"]])
cmd_netMHCpan <- fread(snakemake@input[["cmd_netMHCpan"]])

# Output dirs
dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
dir_DB_out = snakemake@params[["dir_DB_out"]]
dir_IC50 = snakemake@params[["dir_IC50"]]
suppressWarnings(dir.create(dir_DB_out))
suppressWarnings(dir.create(paste0(dir_DB_out, "/plots")))
suppressWarnings(dir.create(paste0(dir_DB_out, "/chunks")))

### ---------------------------- (2) Pre-processing --------------------------------------
# Find all non-empty IC50 prediction outputs
IC50_output <- list.files(paste0(dir_IC50, "/IC50_filtered_peptides"), pattern = ".fst", full.names = T) %>%
  as_tibble() %>%
  mutate(size = file.size(value)) %>%
  filter(size > 40) %>%
  mutate(Peptide_file = str_remove_all(value, ".fst")) %>%
  mutate(Peptide_file = str_split_fixed(Peptide_file, "/IC50_filtered_peptides/", n = 2)[,2]) %>%
  mutate(Affinity_threshold = as.numeric(str_extract(Peptide_file, "[\\d+]+$"))) %>%
  mutate(Peptide_file = str_sub(Peptide_file, start = 1, end = str_length(Peptide_file) - str_length(Affinity_threshold) - 1)) %>%
  mutate(AA = str_split_fixed(Peptide_file, "_", n = 4)[,1]) %>%
  mutate(Nmer = str_split_fixed(Peptide_file, "_", n = 4)[,2]) %>%
  mutate(allele = str_split_fixed(Peptide_file, "_", n = 4)[,3]) %>%
  unite(col = AA_length, AA, Nmer, sep = "_")
IC50_output

# Add non-binder MW_RT filter outputs if necessary
Experiment_design_expanded <- Experiment_design %>%
  as_tibble() %>%
  mutate(`MHC-I_alleles` = ifelse(is.na(`MHC-I_alleles`), "", `MHC-I_alleles`)) %>%
  mutate(Affinity_threshold = ifelse(is.na(Affinity_threshold), "", Affinity_threshold)) %>%
  tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]") %>%
  mutate(`MHC-I_alleles` = str_squish(`MHC-I_alleles`),
         Affinity_threshold = str_squish(Affinity_threshold)) 
Experiment_design_expanded

{
  AA_lengths <- Master_table_expanded %>%
    as_tibble() %>%
    select(N_mers) %>%
    tidyr::expand_grid(AA) %>%
    tidyr::unite(AA_length, AA, N_mers, sep = "_") %>%
    pull(AA_length) %>%
    str_c(collapse = "|")
  
  Mass_list_files <- Experiment_design$Filename %>%
    unique() %>%
    str_c(collapse = "|")
  
  MW_RT_output <- list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched"), pattern = ".fst", full.names = F) %>%
    as_tibble() %>%
    mutate(Mass_list_file = str_extract(value, pattern = Mass_list_files)) %>%
    mutate(enzyme_type = str_split_fixed(value, "_", n = 4)[,1]) %>%
    mutate(AA_length = str_extract(value, AA_lengths)) %>%
    mutate(AA = str_split_fixed(AA_length, "_", n = 2)[,1]) %>%
    mutate(Nmer = str_sub(AA_length, start = 2 + str_length(AA))) %>%
    unite(col = AA_length, AA, Nmer, sep = "_") %>%
    left_join(select(Experiment_design_expanded, Filename, Biological_group, `Output non-binders`), 
              by = c("Mass_list_file" = "Filename")) %>%
    mutate(`Output non-binders` = ifelse(`Output non-binders` == TRUE, TRUE, FALSE)) %>%
    mutate(`Output non-binders` = ifelse(is.na(`Output non-binders`), FALSE, `Output non-binders`)) %>%
    unique()
  MW_RT_output
}

# Select all required peptide mapping tables
{
  Proteomes <- Master_table_expanded$Proteome %>%
    unique() %>%
    str_c(collapse = "|")
  
  pep_map <- list.files(paste0(dir_DB_exhaustive, "/peptide_mapping"), pattern = ".fst", full.names = F) %>%
    as_tibble() %>%
    mutate(enzyme_type = str_split_fixed(value, "_map_", n = 2)[,1]) %>%
    mutate(AA_length = str_extract(value, AA_lengths)) %>%
    mutate(AA = str_split_fixed(AA_length, "_", n = 2)[,1]) %>%
    mutate(Nmer = str_sub(AA_length, start = 2 + str_length(AA))) %>%
    unite(col = AA_length, AA, Nmer, sep = "_") %>%
    mutate(Proteome = str_extract(value, pattern = Proteomes)) %>%
    filter(!str_starts(AA_length, "_")) %>%
    as.data.table() %>%
    split(by = "AA_length")
  pep_map
}

# --------------------------------------- Main --------------------------------------
# x=pep_map[[3]]
DB_reduction_IC50 <- lapply(pep_map, FUN = function(x){
  print(x$AA_length[[1]]) 
  fst::threads_fst(nr_of_threads = Ncpu)
  
  # Load IC50 filtered peptides
  {
    IC50_x <- IC50_output %>%
      filter(AA_length %in% x$AA_length)
    IC50_pep <- IC50_x %>%
      pull(value) %>%
      as.list() %>%
      bettermc::mclapply(mc.cores = Ncpu, mc.progress = F, mc.preschedule = F, FUN = read_fst, as.data.table = TRUE) %>%
      rbindlist()
    
    if (nrow(IC50_pep) == 0) {
      IC50_pep <- data.table(MHC = NA,
                             peptide = NA,
                             'Aff(nM)' = NA,
                             Predicted_binder = F)
    } else {
      setkey(IC50_pep, "peptide")
    }
  }
  
  # Load peptide-protein mappings
  {
    peptide_mapping <- paste0(dir_DB_exhaustive, "/peptide_mapping/", x$value) %>%
      bettermc::mclapply(mc.cores = Ncpu, mc.progress = F,  mc.preschedule = F, FUN = read_fst, as.data.table = TRUE)
    names(peptide_mapping) <- paste(x$enzyme_type, x$Proteome, sep = "|")
    
    peptide_mapping <- rbindlist(peptide_mapping, idcol="file")
    peptide_mapping[, enzyme_type := str_split_fixed(file, "\\|", 2)[,1]]
    peptide_mapping[, Proteome := str_split_fixed(file, "\\|", 2)[,2]]
    peptide_mapping$file <- NULL
    peptide_mapping[, protein := str_split_fixed(protein, "\\.|chunk:", 2)[,1]]
    setkey(peptide_mapping, "peptide")
    peptide_mapping <- unique(peptide_mapping)
    peptide_mapping
  }
  fst::threads_fst(nr_of_threads = Ncpu)
  
  # Load MW RT and define biological groups
  {
    MW_RT_x <- MW_RT_output %>%
      filter(AA_length %in% x$AA_length) %>%
      arrange(Biological_group, `Output non-binders`)
    
    Biological_groups <- unique(MW_RT_x$Biological_group)
    Biological_groups_stats <- vector(mode = "list", length = length(Biological_groups))
  }
  # biol_group=Biological_groups[[2]]
  for (biol_group in Biological_groups) {
    # print(biol_group)
    
    MW_RT_biol_group <- MW_RT_x %>%
      filter(Biological_group %in% biol_group) 
    MW_RT_biol_group_pep <- vector(mode = "list", length = nrow(MW_RT_biol_group))
    
    # Join IC50 with MW_RT information
    for (i in seq_along(MW_RT_biol_group_pep)) {
      # print(i)
      
      # Load 2D filtered peptides
      MW_RT_biol_group_pep[[i]] <- read_fst(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", MW_RT_biol_group$value[[i]]), as.data.table = T)
      setkey(MW_RT_biol_group_pep[[i]], "peptide")
      
      if (MW_RT_biol_group$`Output non-binders`[[i]] == TRUE) {
        MW_RT_biol_group_pep[[i]] <- merge(MW_RT_biol_group_pep[[i]], IC50_pep, all.x = TRUE) %>%
          unique()
        
      } else {
        MW_RT_biol_group_pep[[i]] <- merge(MW_RT_biol_group_pep[[i]], IC50_pep, all.y = TRUE) %>%
          unique() %>%
          lazy_dt() %>%
          filter((!is.na(MW)) & Predicted_binder == TRUE) %>%
          as.data.table()
      }
    }
    ### Report number of predicted binders per mass_list file
    IC50_filtered_stats <- bettermc::mclapply(mc.cores = Ncpu, mc.progress = F,
                                              MW_RT_biol_group_pep, FUN = function(x){
                                                if (nrow(x) > 0 ) {
                                                  out <- x[,length:=str_length(peptide)] %>% 
                                                    split(x$length) %>%
                                                    lapply(function(df){
                                                      df %>%
                                                        filter(!is.na(`Aff(nM)`)) %>%
                                                        filter(Predicted_binder == T) %>%
                                                        select(peptide) %>%
                                                        as_tibble() %>%
                                                        n_distinct()
                                                    }) %>%
                                                    stack() %>%
                                                    rename(IC50_filtered_peptides = values,
                                                           length = ind) %>%
                                                    mutate(length = as.integer(as.character(length)))
                                                } else {
                                                  out <- tibble(IC50_filtered_peptides = 0,
                                                                length = 0)
                                                }
                                              }) %>%
      bind_rows()
    
    # MW_RT_biol_group <- expand_grid(MW_RT_biol_group, IC50_filtered_stats)
    MW_RT_biol_group <- bind_cols(MW_RT_biol_group, IC50_filtered_stats)
    
    # Proceed with non-empty data.tables
    keep <- lapply(lapply(MW_RT_biol_group_pep, dim), `[[`, 1) > 0
    MW_RT_biol_group_pep <- MW_RT_biol_group_pep[keep]
    
    # Save non-empty outputs
    if (length(MW_RT_biol_group_pep) > 0) {
      
      # Join with pep_map
      keep_pep <- rbindlist(MW_RT_biol_group_pep) 
      setkey(keep_pep, peptide) 
      keep_pep <- merge(keep_pep, peptide_mapping, all.x=TRUE) %>%
        unique() %>%
        as_tibble() %>%
        relocate(peptide, enzyme_type, Proteome, protein, MW, RT_pred, MHC, `Aff(nM)`) 
      
      # Save tabular
      fwrite(keep_pep, file = paste0(dir_DB_out, "/", biol_group, ".csv"), sep = ",", append = T, nThread = Ncpu)
      
      # Make fasta headers
      if (minimal_output_headers == T) {
        keep_pep <- keep_pep %>%
          select(peptide) %>%
          unique() %>%
          mutate(header=paste0(biol_group, "_", x$AA_length[[1]], "_", 1:n())) %>%
          as.data.table()
      } else {
        keep_pep <- keep_pep %>%
          unite(col = header, enzyme_type, Proteome, protein, MW, RT_pred, MHC, `Aff(nM)`, remove = T, sep = "|") %>%
          lazy_dt() %>%
          group_by(peptide) %>%
          summarise(header=paste(header,collapse=';')) %>%
          as.data.table()
      }
      # Save fasta
      fasta_chunk <- AAStringSet(x = keep_pep$peptide)
      names(fasta_chunk) <- keep_pep$header
      fasta_chunk <- unique(fasta_chunk)
      writeXStringSet(fasta_chunk, append = T, format = "fasta",
                      filepath = paste0(dir_DB_out, "/", biol_group, ".fasta"))
    }
    # Return stats for a biological group
    Biological_groups_stats[[biol_group]] <- MW_RT_biol_group
  }
  out <- rbindlist(Biological_groups_stats)
  fwrite(out, file = paste0(dir_DB_out, "/chunks/", x$AA_length[[1]], ".csv"))
  return(out)
}) %>%
  bind_rows()

### ---------------------------- (3) Summarize filtering stats --------------------------------------
MW_RT_stats <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status"), full.names = T) %>%
  bettermc::mclapply(mc.cores = Ncpu, mc.progress = F,
                     FUN = fread, nThread = 1, sep = ",") %>%
  rbindlist() %>%
  filter(!is.na(unique_peptides)) %>%
  as_tibble() %>%
  rename(enzyme_type = Splice_type) %>%
  mutate(AA = ifelse(AA == "FALSE" | AA == FALSE, "F", AA)) %>%
  mutate(AA = ifelse(AA == "TRUE" | AA == TRUE, "T", AA)) %>%
  select(-c("value", "Time", "filename", "file", "AA_length", "AA", "PTMs"))  %>%
  unique() 
MW_RT_stats

### Add peptide filtering information from other steps
tmp <- Experiment_design_expanded %>%
  select(Filename, Biological_group, "Output non-binders") %>%
  rename(mass_list = Filename) %>%
  unique()
tmp

Summary_stats <- DB_reduction_IC50 %>%
  as_tibble() %>%
  mutate(length = as.integer(as.character(length))) %>%
  rename(mass_list = Mass_list_file) %>%
  mutate(Proteome = str_extract(value, pattern = str_c(unique(Master_table_expanded$Proteome), 
                                                       collapse = "|"))) %>%
  
  select(-c("value", "AA_length")) 
Summary_stats

Summary_stats <- MW_RT_stats %>%
  left_join(Summary_stats) %>%
  # Fill in NAs that appear due to complete removal of AA_Nmer at 2D filter step
  mutate(IC50_filtered_peptides = ifelse(is.na(IC50_filtered_peptides), 0, IC50_filtered_peptides)) %>%
  select(-c("Biological_group", "Output non-binders")) %>%
  left_join(tmp) %>%
  # Reshape for plotting
  pivot_longer(cols = contains("_peptides"), names_to = c("Filtering_step"), values_to = "# peptides") %>%
  mutate(`log10 # peptides` = log10(`# peptides` + 1)) %>%
  mutate(Filtering_step = str_remove(Filtering_step, "_peptides")) %>%
  mutate(Filtering_step = factor(Filtering_step, levels = c("unique", "mz_matched", "mz_RT_matched", "IC50_filtered"))) %>%
  unique()

### Plot DB size reduction
gg$DB_reduction <- Summary_stats %>%
  # Summarize
  group_by(mass_list, enzyme_type, Proteome, Filtering_step, length) %>%
  summarise(`# peptides` = sum(`# peptides`)) %>%
  mutate(`log10 # peptides` = log10(`# peptides` + 1)) %>%
  mutate(Length = as.factor(length)) %>%
  
  # Plot
  ggplot(aes(y=`log10 # peptides`, x=as.factor(Filtering_step), group=interaction(Proteome, Length))) + 
  geom_line(aes(color=Proteome, linetype = Length), alpha=1) + 
  facet_grid(enzyme_type ~ mass_list) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  guides(color=guide_legend(title="Proteome")) + 
  xlab("Filtering step") + 
  ylab("log10 # peptides") 
gg$DB_reduction

### ---------------------------- (4) Save outputs --------------------------------------
# Plots
for (i in seq_along(gg)) {
  ggsave(gg[[i]], filename = paste0(dir_DB_out, "/plots/", names(gg)[[i]], ".png"), device = "png", 
         width = 20, height = 6, dpi = "retina")
  ggsave(gg[[i]], filename = paste0(dir_DB_out, "/plots/", names(gg)[[i]], ".pdf"), device = "pdf", 
         width = 20, height = 6, dpi = "retina")
}

# Stats
Summary_stats %>%
  # Summarize
  group_by(mass_list, enzyme_type, Proteome, Filtering_step, length) %>%
  summarise(`# peptides` = sum(`# peptides`)) %>%
  mutate(`log10 # peptides` = log10(`# peptides` + 1)) %>%
  pivot_wider(id_cols = c("mass_list", "enzyme_type", "Proteome", "length"), 
              names_from = Filtering_step, values_from = c(`# peptides`)) %>%
  fwrite(sep = ",", append = FALSE, col.names = TRUE,
         file = unlist(snakemake@output[["Summary_stats"]]))
sink()

# Summary_stats %>%
#   group_by(mass_list, enzyme_type, Proteome, Filtering_step, length) %>%
#   summarise(`# peptides` = sum(`# peptides`)) %>%
#   mutate(`log10 # peptides` = log10(`# peptides` + 1)) %>%
#   pivot_wider(id_cols = c("mass_list", "enzyme_type", "Proteome", "length"),
#               names_from = Filtering_step, values_from = c(`# peptides`)) %>%
#   fwrite(sep = ",", append = FALSE, col.names = TRUE,
#          file = paste0(dir_DB_out, "/Summary_stats.csv"))
