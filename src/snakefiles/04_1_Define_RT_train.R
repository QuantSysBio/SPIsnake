### ---------------------------------------------- Define RT train  ----------------------------------------------
# description:
#               
# input:        1. RT calibration datasets: sequence + time (seconds)
#               2. Parameters: n_folds, train_proportion
# output:       
#               A table with a single line per combination of parameters across calibration datasets. 
#               Data split into train/test sets
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))

{
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
  # n_folds = 3
  # proportion_train = 0.9
  # dir_RT_calibration = "data/RT_calibration/"
  # dir_RT_prediction = "results/RT_prediction/"
}
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (1) Read inputs ----------------------------
# Master_table_expanded
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)

# train/test split
dir_RT_calibration = snakemake@params[["dir_RT_calibration"]]
n_folds = as.integer(snakemake@params[["n_folds"]])
proportion_train = as.numeric(snakemake@params[["proportion_train"]])

# Output dirs
dir_RT_prediction = snakemake@params[["dir_RT_prediction"]]
{
  suppressWarnings(
    dir.create(dir_RT_prediction)
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/train"))
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/test"))
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/predict"))
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/AutoRT_models"))
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/plots"))
  )
  suppressWarnings(
    dir.create(paste0(dir_RT_prediction, "/peptide_RT"))
  )
}







### ---------------------------- (2) Split train/test --------------------------------------
# Prepare calibration data
train <- list.files(paste0(dir_RT_calibration), full.names = T) %>%
  lapply(vroom, show_col_types = F) 

names(train) <- list.files(paste0(dir_RT_calibration), full.names = F) %>%
  str_remove_all(".csv")

train <- train %>%
  bind_rows(.id="dataset") %>%
  rename(RT = rt,
         peptide = s) %>%
  select(dataset, peptide, RT) %>%
  split(~dataset, drop = T)

# Save files
for (i in seq_along(train)) {
  for (j in 1:n_folds) {

    train_i = train[[i]] %>%
      ungroup() %>%
      select(peptide, RT) %>%
      unique() %>%
      slice_sample(prop = proportion_train)
    
    test_i = train[[i]] %>%
      filter(!peptide %in% train_i$peptide) %>%
      ungroup() %>%
      select(peptide, RT) %>%
      unique()
    
    ### AutoRT
    # For evaluating prediction error
    train_i %>%
      rename(x=peptide, y=RT) %>%
      select(x, y) %>%
      vroom_write(file = paste0(dir_RT_prediction, "/train/", names(train)[[i]], ".sample", j, ".tsv"), 
                  num_threads = Ncpu, append = F, delim = "\t")
    
    test_i %>%
      rename(x=peptide, y=RT) %>%
      select(x, y) %>%
      vroom_write(file = paste0(dir_RT_prediction, "/test/", names(train)[[i]], ".sample", j, ".tsv"), 
                  num_threads = Ncpu, append = F, delim = "\t")
  }
  # For training on 100% of data
  train[[i]] %>%
    ungroup() %>%
    rename(x=peptide, y=RT) %>%
    select(x, y) %>%
    unique() %>%
    vroom_write(file = paste0(dir_RT_prediction, "/train/", names(train)[[i]], ".tsv"), 
                num_threads = Ncpu, append = F, delim = "\t")
}

### ---------------------------- (3) Define cmds: n-fold cross-validation --------------------------------------
# prefix = ".././workdir/" # for Docker 
prefix = ""
test_AutoRT = paste0(prefix, c(
  list.files("results/RT_prediction/test", full.names = T, pattern = ".tsv")
))
train_AutoRT = paste0(prefix, c(
  # list.files("results/RT_prediction/AutoRT_models_100/", full.names = T, pattern = ".tsv"),
  list.files("results/RT_prediction/train", full.names = T, pattern = ".tsv")
))

out_AutoRT = train_AutoRT %>%
  str_replace(pattern = "results/RT_prediction/train/", replacement = "results/RT_prediction/AutoRT_models/") %>%
  str_remove_all(pattern = fixed(".tsv")) 

# Dirs to store models
out_AutoRT %>%
  str_remove_all(fixed(".././workdir/")) %>%
  lapply(dir.create, showWarnings = FALSE)

# Dirs to store predictions
pred_AutoRT <- out_AutoRT[str_ends(out_AutoRT, ".sample[:digit:]")] %>%
  str_replace(pattern = "results/RT_prediction/AutoRT_models/", replacement = "results/RT_prediction/predict/")

pred_AutoRT %>%
  str_remove_all(fixed(".././workdir/")) %>%
  lapply(dir.create, showWarnings = FALSE)

# Train commands
cmd_AutoRT_train <- c()
for (i in seq_along(train_AutoRT)) {
  cmd_AutoRT_train[i] <- paste("python bin/AutoRT/autort.py train -i", train_AutoRT[i] ,
                               "-o", out_AutoRT[i] ,
                               "-e 40 -b 64 -u m -m bin/AutoRT/models/general_base_model/model.json",
                               "--add_ReduceLROnPlateau --early_stop_patience 10")
}

# Predict commands
cmd_AutoRT_test <- c()
for (i in seq_along(pred_AutoRT)) {
  cmd_AutoRT_test[i] <- paste("python bin/AutoRT/autort.py predict --test", test_AutoRT[i],
                              "-s", paste0(out_AutoRT[i], "/model.json"),
                              "-o", pred_AutoRT[i])
}

dataset_name <- str_split_fixed(train_AutoRT, pattern = "results/RT_prediction/", 2)[,2] %>%
  str_split_fixed(pattern = "/", 2)
dataset_name <- dataset_name[,2] %>%
  str_remove_all("/") %>%
  str_remove_all(".tsv")

### ---------------------------- (4) Export cmds --------------------------------------
{
  # data.frame(RT_dataset = dataset_name,
  #            cmd = cmd_AutoRT_train)  %>%
  #   vroom_write(delim = ",", append = FALSE,
  #               file = "results/RT_prediction/cmd_AutoRT_train.csv")
  # 
  # data.frame(number = 1:length(cmd_AutoRT_test),
  #            cmd = cmd_AutoRT_test) %>%
  #   vroom_write(delim = ",", append = FALSE,
  #               file = "results/RT_prediction/cmd_AutoRT_test.csv")
}
data.frame(RT_dataset = dataset_name,
           cmd = cmd_AutoRT_train)  %>%
  vroom_write(delim = ",", append = FALSE,
              file = unlist(snakemake@output[["cmd_AutoRT_train"]]))

data.frame(number = 1:length(cmd_AutoRT_test),
           cmd = cmd_AutoRT_test) %>%
  vroom_write(delim = ",", append = FALSE,
              file =  unlist(snakemake@output[["cmd_AutoRT_test"]]))

