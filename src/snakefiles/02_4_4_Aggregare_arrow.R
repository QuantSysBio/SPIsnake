### ---------------------------------------------- Agrregate arrow ----------------------------------------------
# description:  Create unique arrow outputs across chunks
#               
# input:        1. Generate_peptides rules are done
#               2. Peptide sequences are saved as arrow datasets
# output:       
#               - Unique peptides as arrow datasets collected across chunks
#               
# author:       Yehor Horokhovskyi

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, split = TRUE)
}
cat(as.character(Sys.time()), " - ", "Started R", "\n")
cat(as.character(Sys.time()), " - ", R.utils::getHostname.System(), "\n")

suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(bettermc))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read input file and extract info ----------------------------
if (exists("snakemake")) {
  arrow_batch_definition <- fread(snakemake@input[["arrow_batch_definition"]])
  glimpse(arrow_batch_definition)
  # DB_PTM_mz <- fread(snakemake@input[["DB_PTM_mz"]])
  
  # CPUs
  active_Slurm <- ifelse(system("echo $SLURM_JOB_ID") != 0, TRUE, FALSE)
  print(active_Slurm)
  Ncpu <- snakemake@params[["max_cpus"]]
  print(paste0("Ncpu: ", Ncpu))
  cl <- parallel::makePSOCKcluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  
  # Folders
  dir_DB_PTM_mz = snakemake@params[["dir_DB_PTM_mz"]]
  dir_DB_out = snakemake@params[["dir_DB_out"]]
  duckdb_RAM = snakemake@params[["duckdb_RAM"]]
  duckdb_max_retries = snakemake@params[["duckdb_max_retries"]]
  duckdb_max_filesize = snakemake@params[["duckdb_max_filesize"]]
  
  # Wildcard
  filename <- snakemake@output[[1]]
} else {
  arrow_batch_definition <- fread("results/DB_out/arrow_batch_definition.csv")
  # DB_PTM_mz <- fread("results/DB_out/DB_PTM_mz.csv")
  
  # CPUs
  Ncpu <- min(parallelly::availableCores() - 1)
  print(paste0("number of CPUs: ", Ncpu))
  cl <- parallel::makePSOCKcluster(Ncpu)
  cl <- parallelly::autoStopCluster(cl, debug=T)
  
  ### duckdb settings
  duckdb_RAM <- 0.8
  duckdb_max_retries = 5
  duckdb_max_filesize = 2
  
  # Folders
  dir_DB_PTM_mz <- "results/DB_PTM_mz/" 
  dir_DB_out = "results/DB_out/" 
  
  # Wildcard
  filename <- "results/DB_out/.peptides_13.done"
}

##№ duckdb parameters
duckdb_max_filesize <- duckdb_max_filesize * (2^30)
cat(as.character(Sys.time()), " - ", "duckdb max filesize: ", duckdb_max_filesize, "\n")

aggregation_batch_i <- str_remove(filename, ".done")
aggregation_batch_i <- str_split_fixed(aggregation_batch_i, "results/DB_out/.peptides_", 2)[,2][[1]]
aggregation_batch_i

duckdb_temp_dir <- paste0(dir_DB_out, "/duckdb/tmp/unique_peptides_", aggregation_batch_i, "/")
table_name = paste0(dir_DB_out, "/duckdb/databases/", "duckdb_aggregate_unique_", aggregation_batch_i)
Max_RAM <- 240
timeout <- 90
duckdb_RAM <- Max_RAM * duckdb_RAM
duckdb_max_retries = 20

### ---------------------------- (2) Aggregate arrow across chunks --------------------------------------
dir.create(paste0(dir_DB_out, "/arrow/"), recursive = T, showWarnings = F)
dir.create(paste0(dir_DB_out, "/duckdb/databases/"), recursive = T, showWarnings = F)

if (dir.exists(duckdb_temp_dir)) {
  unlink(duckdb_temp_dir, recursive = T, force = T)
} else {
  dir.create(duckdb_temp_dir, recursive = T, showWarnings = T)
}

### ------------------------------- dbplyr solution (duckdb) -------------------------------
### Check no coersion to NA
arrow_batch_definition[chunk == "index=NA", index := "NA"]

### how many CPUs to use
aggregation_CPU <- ifelse(nchar(arrow_batch_definition$index[[1]]) == 1, Ncpu, ceiling(Ncpu / 10))
cat(as.character(Sys.time()), " - ", "CPUs to use for aggregation:", aggregation_CPU, "\n")

### Proteomes 
proteomes <- fread("Master_table.csv") %>%
  as_tibble() %>%
  pull(Proteome)

### Which chunks to process
DB_groups <- arrow_batch_definition[aggregation_batch == aggregation_batch_i,]
glimpse(DB_groups)

chunks <- nrow(DB_groups)
cat(as.character(Sys.time()), " - ", "Defined chunks for processing:", chunks, "\n")
done <- bettermc::mclapply(X = 1:chunks, 
                           mc.preschedule = F, mc.cores = aggregation_CPU, 
                           mc.cleanup = T, mc.retry = 3, 
                           FUN = function(i){
                             cat(as.character(Sys.time()), " - ", i, "/", chunks,"\n")
                             DB_groups_i <- DB_groups[i,] %>%
                               as_tibble()
                             
                             # Let each thread use a new connection
                             table_name <- paste0(table_name, "_", i)
                             
                             ### Connect database
                             cat(as.character(Sys.time()), " - ", "Start creating duckdb", "\n")
                             conn <- DBI::dbConnect(
                               duckdb(), 
                               dbname = paste0("arrow_unique.duckdb.", aggregation_batch_i, "_", i),
                               config=list("memory_limit"= paste0(duckdb_RAM, "GB"),
                                           "temp_directory" = paste0(duckdb_temp_dir, "_", i)))
                             cat(as.character(Sys.time()), " - ", "Done creating duckdb", "\n")
                             
                             DB_duck <- get_distict_duckdb(input_path = paste0(dir_DB_PTM_mz, "peptide_seqences/"), 
                                                           DB_groups_i = DB_groups_i,
                                                           table_name = table_name,
                                                           conn = conn)
                             for (p in proteomes) {
                               DB_duck %>%
                                 filter(proteome == local(p)) %>%
                                 to_arrow() %>%
                                 group_by(index, length, proteome, enzyme, MiSl) %>%
                                 write_dataset(path = paste0(dir_DB_out, "/arrow/"),
                                               existing_data_behavior = "overwrite",
                                               format = "parquet",
                                               max_partitions = 10240L,
                                               max_rows_per_file = as.integer(5 * 10^7),
                                               compression = "lz4")
                             }
                             rm(DB_duck)
                             cat(as.character(Sys.time()), " - ", "Done save arrow", "\n")
                             return(T)
                           }) %>%
  suppressMessages()
cat(as.character(Sys.time()), " - ", "Done saving aggregated arrow dataset", "\n")

# Remove duckdb temp dir
if (dir.exists(duckdb_temp_dir)) {
  unlink(duckdb_temp_dir, recursive = T, force = T)
}
file.remove(paste0("arrow_unique.duckdb.", aggregation_batch_i))

### Log
if (exists("snakemake")) {
  SPIsnake_log()
  sink()
}
