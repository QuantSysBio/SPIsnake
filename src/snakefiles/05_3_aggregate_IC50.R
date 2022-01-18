### ---------------------------------------------- Aggregate IC50  ----------------------------------------------
# description:
#               
# input:        1. metadata: Master_table_expanded, cmd_netMHCpan
#               2. Candidate peptide binders filtered from netMHCpan output
#               
# output:       FASTA file per Biological group
#               multimappers 
#                
#               
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))
print(sessionInfo())

# {
#   ### setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
#   Master_table_expanded <- vroom("results/DB_exhaustive/Master_table_expanded.csv", show_col_types = FALSE)
#   Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE) 
#   dir_DB_exhaustive = "results/DB_exhaustive/"
#   dir_DB_PTM_mz = "results/DB_PTM_mz"
#   dir_IC50 = "results/IC50/"
#   dir_DB_out = "results/DB_out"
#   cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
# }

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
gg <- list()

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

# create temporary directory for vroom
{
  Sys.getenv("TMPDIR") %>% print()
  Sys.getenv("VROOM_TEMP_PATH") %>% print()
  
  vroom_dir = "/tmp/vroom"
  suppressWarnings(dir.create(vroom_dir))
  Sys.setenv(VROOM_TEMP_PATH = vroom_dir)
  Sys.getenv("VROOM_TEMP_PATH") %>%print()
  
  tmp_file = tempfile()
  print(tmp_file)
}

### ---------------------------- (1) Read inputs ----------------------------
# Experiment_design
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE) 
# %>%  mutate(`Output non-binders` = ifelse(Filename %in% c("MeV_BLCL_allFractions", "MeV_MA0009-BE08_allFractions"), T, F))
Master_table_expanded <- vroom(snakemake@input[["Master_table_expanded"]], show_col_types = FALSE)
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)

# Output dirs
dir_DB_exhaustive = snakemake@params[["dir_DB_exhaustive"]]
dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
dir_DB_out = snakemake@params[["dir_DB_out"]]
dir_IC50 = snakemake@params[["dir_IC50"]]
suppressWarnings(dir.create(dir_DB_out))
suppressWarnings(dir.create(paste0(dir_DB_out, "/plots")))

### ---------------------------- (2) Pre-processing --------------------------------------
# Find all non-empty IC50 prediction outputs
IC50_output <- list.files(paste0(dir_IC50, "/IC50_filtered_peptides"), pattern = ".csv.gz", full.names = T) %>%
  as_tibble() %>%
  mutate(size = file.size(value)) %>%
  filter(size > 40) %>%
  mutate(Peptide_file = str_remove_all(value, ".csv.gz")) %>%
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
  tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]") %>%
  mutate(`MHC-I_alleles` = str_squish(`MHC-I_alleles`),
         Affinity_threshold = str_squish(Affinity_threshold)) 
Experiment_design_expanded

# Remove the if condition?
# if (TRUE %in% Experiment_design_expanded$`Output non-binders`) 
{
  Mass_list_files <- Experiment_design$Filename %>%
    unique() %>%
    str_c(collapse = "|")
  
  MW_RT_output <- list.files(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched"), pattern = ".csv.gz", full.names = F) %>%
    as_tibble() %>%
    mutate(Mass_list_file = str_extract(value, pattern = Mass_list_files)) %>%
    mutate(enzyme_type = str_split_fixed(value, "_", n = 4)[,1]) %>%
    mutate(AA = str_split_fixed(value, "_", n = 4)[,2]) %>%
    mutate(Nmer = str_split_fixed(value, "_", n = 4)[,3]) %>%
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
  
  pep_map <- list.files(paste0(dir_DB_exhaustive, "/peptide_mapping"), pattern = ".csv.gz", full.names = F) %>%
    as_tibble() %>%
    mutate(enzyme_type = str_split_fixed(value, "_map_", n = 2)[,1]) %>%
    mutate(AA_length = str_split_fixed(value, "_map_", n = 2)[,2]) %>%
    mutate(AA = str_split_fixed(AA_length, "_", n = 2)[,1]) %>%
    mutate(Nmer = str_split_fixed(AA_length, "_", n = 3)[,2]) %>%
    unite(col = AA_length, AA, Nmer, sep = "_") %>%
    mutate(Proteome = str_extract(value, pattern = Proteomes)) %>%
    filter(!str_starts(AA_length, "_")) %>%
    # filter(AA_length %in% c(MW_RT_output$AA_length[MW_RT_output$`Output non-binders` == T], IC50_output$AA_length)) %>%
    as.data.table() %>%
    split(by = "AA_length")
  pep_map
}

# ----------------------------- DEV: ON ----------------------------
DB_reduction_IC50 <- lapply(pep_map, FUN = function(x){
  print(x$AA_length[[1]])
  
  # Load IC50 filtered peptides
  {
    IC50_x <- IC50_output %>%
      filter(AA_length %in% x$AA_length)
    IC50_pep <- IC50_x %>%
      pull(value) %>%
      vroom(show_col_types = F, num_threads = Ncpu, delim = ",") %>%
      as.data.table()
    
    if (nrow(IC50_pep) == 0) {
      IC50_pep <- data.table(MHC = NA,
                             peptide = NA,
                             'Aff(nM)' = NA)
      }
  }
  
  # Load peptide-protein mappings
  {
    peptide_mapping <- paste0(dir_DB_exhaustive, "/peptide_mapping/", x$value) %>%
      bettermc::mclapply(mc.cores = 1, mc.progress = F,
                         FUN = vroom, show_col_types = F, num_threads = 1, delim = ",")
    names(peptide_mapping) <- paste(x$enzyme_type, x$Proteome, sep = "|")
    
    peptide_mapping <- rbindlist(peptide_mapping, idcol="file")
    # peptide_mapping[, c("enzyme_type", "Proteome") := tstrsplit(file, "|", fixed=TRUE)]
    peptide_mapping[, enzyme_type := str_split_fixed(file, "\\|", 2)[,1]]
    peptide_mapping[, Proteome := str_split_fixed(file, "\\|", 2)[,2]]
    peptide_mapping$file <- NULL
    peptide_mapping[, protein := str_split_fixed(protein, "\\|chunk:", 2)[,1]]
    peptide_mapping <- unique(peptide_mapping)
    peptide_mapping
  }
  
  # Load MW RT and define biological groups
  {
    MW_RT_x <- MW_RT_output %>%
      filter(AA_length %in% x$AA_length) %>%
      arrange(Biological_group, `Output non-binders`)
    
    Biological_groups <- unique(MW_RT_x$Biological_group)
    Biological_groups_stats <- vector(mode = "list", length = length(Biological_groups))
  }
  for (biol_group in Biological_groups) {
    
    MW_RT_biol_group <- MW_RT_x %>%
      filter(Biological_group %in% biol_group) %>%
      mutate(IC50_filtered_peptides = NA)
    MW_RT_biol_group_pep <- vector(mode = "list", length = nrow(MW_RT_biol_group))
    
    # Join IC50 with MW_RT information
    for (i in seq_along(MW_RT_biol_group_pep)) {
      if (MW_RT_biol_group$`Output non-binders`[[i]] == TRUE) {
        MW_RT_biol_group_pep[[i]] <- vroom(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", MW_RT_biol_group$value[[i]]), 
                                           show_col_types = F, num_threads = Ncpu, delim = ",") %>%
          lazy_dt() %>%
          left_join(IC50_pep) %>%
          as.data.table()
      } else {
        MW_RT_biol_group_pep[[i]] <- IC50_pep %>%
          lazy_dt() %>%
          left_join(y = vroom(paste0(dir_DB_PTM_mz, "/unique_peptides_mz_RT_matched/", MW_RT_biol_group$value[[i]]), 
                              show_col_types = F, num_threads = Ncpu, delim = ",")) %>%
          filter(!is.na(MW)) %>%
          unique() %>%
          as.data.table()
      }
    }
    # Report number of predicted binders per mass_list file
    MW_RT_biol_group$IC50_filtered_peptides <- bettermc::mclapply(mc.cores = Ncpu, mc.progress = F,
                                                                         MW_RT_biol_group_pep, FUN = function(x){
      x %>%
        filter(!is.na(`Aff(nM)`)) %>%
        select(peptide) %>%
        as_tibble() %>%
        n_distinct()
    }) %>%
      unlist()
    
    # Proceed with non-empty data.tables
    keep <- lapply(lapply(MW_RT_biol_group_pep, dim), `[[`, 1) > 0
    MW_RT_biol_group_pep <- MW_RT_biol_group_pep[keep]
    
    # Save non-empty outputs
    if (length(MW_RT_biol_group_pep) > 0) {
      
      # Join with pep_map
      keep_pep <- MW_RT_biol_group_pep %>%
        bind_rows() %>%
        unique() %>%
        as_tibble() %>%
        left_join(peptide_mapping) %>%
        relocate(peptide, enzyme_type, Proteome, protein, MW, RT_pred, MHC, `Aff(nM)`) %>%
        as.data.table() 
      
      # Save tabular
      vroom_write(keep_pep, file = paste0(dir_DB_out, "/", biol_group, ".csv"), delim = ",", append = T, num_threads = Ncpu)
      
      # Make fasta headers
      keep_pep <- keep_pep %>%
        unite(col = header, enzyme_type, Proteome, protein, MW, RT_pred, MHC, `Aff(nM)`, remove = T, sep = "|") %>%
        lazy_dt() %>%
        group_by(peptide) %>%
        summarise(header=paste(header,collapse=';')) %>%
        as.data.table()
      keep_pep
      
      # Save fasta
      fasta_chunk <- AAStringSet(x = keep_pep$peptide)
      names(fasta_chunk) <- keep_pep$header
      writeXStringSet(fasta_chunk, append = T, format = "fasta",
                      filepath = paste0(dir_DB_out, "/", biol_group, ".fasta"))
    }
    # Return stats for a biological group
    Biological_groups_stats[[biol_group]] <- MW_RT_biol_group
  }
  return(rbindlist(Biological_groups_stats))
}) %>%
  bind_rows()


### ---------------------------- (3) Summarize filtering stats --------------------------------------
MW_RT_stats <- list.files(paste0(dir_DB_PTM_mz, "/chunk_aggregation_status"), full.names = T) %>%
  bettermc::mclapply(mc.cores = Ncpu, mc.progress = F,
                     FUN = vroom, show_col_types = F, num_threads = 1, delim = ",") %>%
  rbindlist() %>%
  select(-c("value", "file", "PTMs","Time", "mass_list_file")) %>%
  filter(!is.na(unique_peptides)) %>%
  unique() %>%
  as_tibble() %>%
  rename(enzyme_type = Splice_type)
MW_RT_stats

### Add peptide filtering information from other steps
Summary_stats <- DB_reduction_IC50 %>%
  as_tibble() %>%
  rename(mass_list = Mass_list_file) %>%
  mutate(Proteome = str_extract(value, pattern = str_c(unique(Master_table_expanded$Proteome), 
                                                       collapse = "|"))) %>%
  full_join(MW_RT_stats) %>%
  pivot_longer(cols = contains("_peptides"), names_to = c("Filtering_step"), values_to = "# peptides") %>%
  mutate(`log10 # peptides` = log10(`# peptides` + 1),
         Length = str_split_fixed(AA_length, "_", 2)[,2]) %>%
  mutate(Filtering_step = str_remove(Filtering_step, "_peptides")) %>%
  mutate(Filtering_step = factor(Filtering_step, levels = c("unique", "mz_matched", "mz_RT_matched", "IC50_filtered")))

### Plot DB size reduction
gg$DB_reduction <- Summary_stats %>%
  # Summarize
  group_by(mass_list, enzyme_type, Proteome, Filtering_step, Length) %>%
  summarise(`# peptides` = sum(`# peptides`)) %>%
  mutate(`log10 # peptides` = log10(`# peptides` + 1)) %>%
  mutate(Length = as.factor(Length)) %>%
  
  # Plot
  ggplot(aes(y=`log10 # peptides`, x=as.factor(Filtering_step), group=interaction(Proteome, Length))) + 
  geom_line(aes(color=Proteome, linetype = Length), alpha=1) + 
  facet_grid(enzyme_type ~ mass_list) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=12, hjust=0.5, vjust=0.5)) + 
  guides(color=guide_legend(title="Proteome")) + 
  xlab("Filtering step") + 
  ylab("log10 # peptides") 
gg$DB_reduction

### ---------------------------- (4) Save outputs --------------------------------------
# Plots
for (i in seq_along(gg)) {
  ggsave(gg[[i]], filename = paste0(dir_DB_out, "/plots/", names(gg)[[i]], ".png"), device = "png", 
         width = 20, height = 6, dpi = "retina")
}

# Stats
Summary_stats %>%
  vroom_write(delim = ",", append = FALSE, col_names = TRUE,
              file = unlist(snakemake@output[["Summary_stats"]]))