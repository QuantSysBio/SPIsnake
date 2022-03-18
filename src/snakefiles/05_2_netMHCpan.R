### ---------------------------------------------- Predict_MHC-I_affinity ----------------------------------------------
# description:  Predict MHC-I affinity for a given combination of alleles and save the filtered peptides 
#               
# input:        1. Peptide sequences generated in chunks
#               2. Parameters: Master_table_expanded
# output:       
#               Output files are saved in chunks that depend on first N=index_length letters of a peptide
#               - Unique sequences per experiment .csv.gz
#               - Unique peptides MW .csv.gz
#               - Unique peptides after m/z and RT matching per experiment .csv.gz
#               - Peptide-mass_list matching
#               - Peptide filtering stats .csv
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log, split = TRUE)

suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fst))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))

source("src/snakefiles/functions.R")
print(sessionInfo())

### ---------------------------- (1) Read input file and extract info ----------------------------
# {
#   ### setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
#   dir_IC50 = "results/IC50/"
#   cmd_netMHCpan <- vroom(paste0(dir_IC50, "cmd_netMHCpan.csv"), show_col_types = FALSE)
#   Experiment_design <- vroom("data/Experiment_design.csv", delim = ",")
#   cmd_block_i = "results/IC50/Seq_stats/1.csv"
#   fst_compression = 100
# }

# Experiment_design
cmd_netMHCpan <- vroom(snakemake@input[["cmd_netMHCpan"]], show_col_types = FALSE)
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)

# Wildcard
cmd_block_i <- snakemake@output[[1]]

# Output dirs
dir_IC50 = snakemake@params[["dir_IC50"]]
suppressWarnings(dir.create(paste0(dir_IC50, "/Seq_stats")))
suppressWarnings(dir.create(paste0(dir_IC50, "/IC50_filtered_peptides")))

# CPU
Ncpu = availableCores()
setDTthreads(Ncpu)

# fst compression level
fst_compression = as.integer(snakemake@params[["fst_compression"]])

### ---------------------------------------------- Predict binding --------------------------------------------------
print(head(Experiment_design))

IC50_aggregation_table <- Experiment_design %>%
  filter(!(is.na(`MHC-I_alleles`) | is.na(Affinity_threshold))) %>%
  tidyr::separate_rows(`MHC-I_alleles`, Affinity_threshold, sep = "[|]") %>%
  mutate(allele = str_squish(`MHC-I_alleles`),
         Affinity_threshold = str_squish(Affinity_threshold)) %>%
  select(Filename, allele, Affinity_threshold) 

# Select commands
cmd_block_i <- str_split_fixed(cmd_block_i, "/IC50/Seq_stats/", 2)[,2] %>%
  str_remove(".csv")
cmd_netMHCpan <- cmd_netMHCpan %>%
  mutate(output_file = str_split_fixed(cmds, pattern = " > ", 2)[,2]) %>%
  mutate(output_file = str_remove_all(output_file, " -v")) %>%
  left_join(IC50_aggregation_table)

# Deselect commands that have already been executed
ready <- list.files(paste0(dir_IC50, "netMHCpan_output"), full.names = T) %>%
  as_tibble() %>%
  mutate(file_size = file.size(value)) %>%
  filter(file_size > 5000) %>%
  mutate(value = str_split_fixed(value, "/netMHCpan_output/", 2)[,2]) %>%
  pull(value) %>%
  str_c(collapse = "|")

cmd_netMHCpan_todo <- cmd_netMHCpan %>% 
  filter(cmd_block == cmd_block_i)

if (!ready == "") {
  cmd_netMHCpan_todo <- cmd_netMHCpan_todo %>%
    filter(!str_ends(output_file, ready))
}
cmd_netMHCpan_todo[1,] %>% t()

# Run netMHCpan-4.1
cmds <- cmd_netMHCpan_todo %>%
  pull(cmds) %>%
  unique() %>%
  as.list()

bettermc::mclapply(cmds, mc.cores = Ncpu, mc.cleanup=T, mc.preschedule=F, function(x){
  print(x)
  system(x, intern = T)
})
print("Predicted MHC affinity")

### ---------------------------------------------- Select MHC-I binders --------------------------------------------------
output_files <- cmd_netMHCpan %>% 
  filter(cmd_block == cmd_block_i) %>%
  select(output_file, Affinity_threshold) %>%
  unique()

output_files_unique <- unique(output_files$output_file)
  
IC50_filter_stats <- bettermc::mclapply(X = output_files_unique, mc.cores = Ncpu, mc.retry = 3, mc.cleanup=T, mc.preschedule=F, FUN = function(x){
  
  if (file.exists(x)) {
    # Parameters
    Affinity_threshold <- cmd_netMHCpan$Affinity_threshold[cmd_netMHCpan$output_file == x] %>%
      unique()
    IC50_filter_stats_AT <- vector(mode = "list", length = length(Affinity_threshold))
    
    for (AT in Affinity_threshold) {
      
      # Read and filter netMHCpan output
      binders_df <- vroom_lines(file = x, 
                                num_threads = 1,
                                skip_empty_rows = TRUE,
                                skip = 55)
      if (length(binders_df) > 0) {
        binders_df <- binders_df %>%
          str_replace_all(pattern = "[[:space:]]+", " ") %>%
          str_split_fixed(pattern = " ", n = Inf) %>%
          as.data.table() %>%
          lazy_dt() %>%
          rename(MHC = V3, peptide = V4, `Aff(nM)` = V17) %>%
          select(MHC, peptide, `Aff(nM)`) %>%
          mutate(`Aff(nM)` = as.numeric(`Aff(nM)`)) %>%
          filter((!is.na(`Aff(nM)`))) %>%
          mutate(Predicted_binder = ifelse(`Aff(nM)` <= as.numeric(AT), TRUE, FALSE)) %>%
          mutate(peptide = str_remove_all(peptide, "X")) %>%
          as.data.table()
        print(head(binders_df))
      } else {
        binders_df = data.table(MHC = character(),
                                peptide = character(),
                                `Aff(nM)` = numeric(),
                                Predicted_binder = logical())
        print("Warning: empty input!")
        print("Make sure that allele is present in netMHCpan-4.1/data/allelenames")
      }
      output_name <- x %>%
        str_remove(".tsv.txt") %>%
        str_replace(pattern = "netMHCpan_output", replacement = "IC50_filtered_peptides")
      
      # Save filtered peptides
      write_fst(binders_df, path = paste0(output_name, "_", AT,".fst"), compress = fst_compression)
      
      # Save stats
      IC50_filter_stats_AT[[AT]] <- cmd_netMHCpan %>%
        filter(output_file == x & Affinity_threshold == AT) %>%
        select(- Filename) %>%
        mutate(IC50_filtered_peptides = length(which(binders_df$Predicted_binder == T)))
    }
    IC50_filter_stats_AT <- bind_rows(IC50_filter_stats_AT)
  } else {
    IC50_filter_stats_AT <- tibble(file = character(),
                                   size = double(),
                                   N_mer = double(),
                                   allele = character(),
                                   cmd_block = double(),
                                   cmds = character(),
                                   output_file = character(),
                                   Affinity_threshold = character(),
                                   IC50_filtered_peptides = double())
  }
  return(IC50_filter_stats_AT)
}) %>%
  bind_rows()
  
### ---------------------------------------------- (4) Export --------------------------------------------------
# IC50_filter_stats %>%
#   vroom_write(append = FALSE, col_names = TRUE, delim = ",",
#               file = "results/IC50/Seq_stats/1.csv")

IC50_filter_stats %>%
  vroom_write(append = FALSE, col_names = TRUE, delim = ",",
              file = unlist(snakemake@output[["IC50_filter_stats"]]))
sink()
