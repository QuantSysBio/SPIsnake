### ---------------------------------------------- Define RT train  ----------------------------------------------
# description:  Evaluate retention time (RT) prediction accuracy.
#               If achrom algorithm is used, predictions are made at once
#               If AutoRT algorithm is used, commands are defined
#               
# input:        1. RT calibration datasets: sequence + time (seconds)
#               2. Parameters: n_folds, train_proportion
# output:       
#               - A table with a single line per combination of parameters across calibration datasets. 
#               - Data split into train/test sets
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
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(vroom))

# {
#   ### setwd("/home/yhorokh/SNAKEMAKE/QSB_SPIsnake/")
#   Experiment_design <- vroom("data/Experiment_design.csv", show_col_types = FALSE)
#   method = "achrom"
#   n_folds = 10
#   proportion_train = 0.9
#   dir_RT_calibration = "data/RT_calibration/"
#   dir_RT_prediction = "results/RT_prediction/"
# }

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### CPUs
Ncpu = availableCores()
cl <- parallel::makeForkCluster(Ncpu)
cl <- parallelly::autoStopCluster(cl)
setDTthreads(Ncpu)

### ---------------------------- (1) Read inputs ----------------------------
# Experiment_design
Experiment_design <- vroom(snakemake@input[["Experiment_design"]], show_col_types = FALSE)

# Retention time prediction method
method = as.character(snakemake@params[["method"]])

# train/test split
n_folds = as.integer(snakemake@params[["n_folds"]])
proportion_train = as.numeric(snakemake@params[["proportion_train"]])
dir_RT_calibration = snakemake@params[["dir_RT_calibration"]]

# Output dirs
dir_RT_prediction = snakemake@params[["dir_RT_prediction"]]
{
  suppressWarnings(dir.create(dir_RT_prediction))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/train")))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/test")))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/predict")))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/RT_models")))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/plots")))
  suppressWarnings(dir.create(paste0(dir_RT_prediction, "/peptide_RT")))
}

### ---------------------------- (2) Split train/test --------------------------------------
# Peptide combinations
AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

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
    
    AA_train <- expand_grid(AA, 
                            peptide = train[[i]]$peptide) %>%
      mutate(detected = str_detect(peptide, AA)) %>%
      select(AA, detected) %>%
      as.data.frame() %>%
      table() %>%
      as.data.frame() %>%
      as_tibble() %>%
      filter(detected == T) %>%
      filter(Freq > 0)
    
    {
      # Make sure that every dataset's AA is in the training data (achrom requirement)
      AA_presence <- expand_grid(AA, 
                                 peptide = train_i$peptide) %>%
        mutate(detected = str_detect(peptide, AA)) %>%
        select(AA, detected) %>%
        as.data.frame() %>%
        table() %>%
        as.data.frame() %>%
        as_tibble() %>%
        filter(detected == TRUE) %>%
        filter(AA %in% AA_train$AA)
      
      while (0 %in% AA_presence$Freq) {
        
        train_i = train[[i]] %>%
          ungroup() %>%
          select(peptide, RT) %>%
          unique() %>%
          slice_sample(prop = proportion_train)
        
        AA_presence <- expand_grid(AA, 
                                   peptide = train_i$peptide) %>%
          mutate(detected = str_detect(peptide, AA)) %>%
          select(AA, detected) %>%
          as.data.frame() %>%
          table() %>%
          as.data.frame() %>%
          as_tibble() %>%
          filter(detected == TRUE) %>%
          filter(AA %in% AA_train$AA)
      }
      }
    
    test_i = train[[i]] %>%
      filter(!peptide %in% train_i$peptide) %>%
      ungroup() %>%
      select(peptide, RT) %>%
      unique()
    
    ### For evaluating prediction error
    # AutoRT
    if (method == "AutoRT") {
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
      
      # For training on 100% of data
      train[[i]] %>%
        ungroup() %>%
        rename(x=peptide, y=RT) %>%
        select(x, y) %>%
        unique() %>%
        vroom_write(file = paste0(dir_RT_prediction, "/train/", names(train)[[i]], ".tsv"), 
                    num_threads = Ncpu, append = F, delim = "\t")
    }
    # Achrom
    if (method == "achrom") {
      vroom_write(train_i, file = paste0(dir_RT_prediction, "/train/", names(train)[[i]], ".sample", j, ".tsv"), 
                  num_threads = Ncpu, append = F)
      vroom_write(test_i, file = paste0(dir_RT_prediction, "/test/", names(train)[[i]], ".sample", j, ".tsv"), 
                  num_threads = Ncpu, append = F)
      # For training on 100% of data
      train[[i]] %>%
        ungroup() %>%
        select(peptide, RT) %>%
        unique() %>%
        vroom_write(file = paste0(dir_RT_prediction, "/train/", names(train)[[i]], ".tsv"), 
                    num_threads = Ncpu, append = F, delim = "\t")
    }
  }
}

### ---------------------------- (3) Define cmds: n-fold cross-validation --------------------------------------
# prefix = ".././workdir/" # for Docker 
prefix = ""
test_RT = paste0(prefix, c(
  list.files("results/RT_prediction/test", full.names = T, pattern = ".tsv")
))
train_RT = paste0(prefix, c(
  # list.files("results/RT_prediction/RT_models_100/", full.names = T, pattern = ".tsv"),
  list.files("results/RT_prediction/train", full.names = T, pattern = ".tsv")
))

out_RT = train_RT %>%
  str_replace(pattern = "results/RT_prediction/train/", replacement = "results/RT_prediction/RT_models/") %>%
  str_remove_all(pattern = fixed(".tsv")) 

# Dirs to store models
out_RT %>%
  str_remove_all(fixed(".././workdir/")) %>%
  lapply(dir.create, showWarnings = FALSE)

# Dirs to store predictions
pred_RT <- out_RT[str_ends(out_RT, ".sample[:digit:]*")] %>%
  str_replace(pattern = "results/RT_prediction/RT_models/", replacement = "results/RT_prediction/predict/")

pred_RT %>%
  str_remove_all(fixed(".././workdir/")) %>%
  lapply(dir.create, showWarnings = FALSE)

# AutoRT
if (method == "AutoRT") {
  # Train commands
  cmd_RT_train <- c()
  for (i in seq_along(train_RT)) {
    cmd_RT_train[i] <- paste("python /bin/AutoRT/autort.py train -i", train_RT[i] ,
                             "-o", out_RT[i] ,
                             "-e 40 -b 64 -u m -m /bin/AutoRT/models/general_base_model/model.json",
                             "--add_ReduceLROnPlateau --early_stop_patience 10")
  }
  
  # Predict commands
  cmd_RT_test <- c()
  for (i in seq_along(pred_RT)) {
    cmd_RT_test[i] <- paste("python /bin/AutoRT/autort.py predict --test", 
                            paste0(str_replace_all(pred_RT[i], pattern = "results/RT_prediction/predict/", replacement = "results/RT_prediction/test/"), ".tsv"),
                            "-s", 
                            paste0(str_replace_all(pred_RT[i], pattern = "results/RT_prediction/predict/", replacement = "results/RT_prediction/RT_models/"), "/model.json"),
                            "-o", pred_RT[i])
  }
} else if (method == "achrom") {
  
  use_condaenv("R_env_reticulate")
  pyteomics <- import("pyteomics")
  
  py_run_string("
import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None
")
  
  # Train
  for (i in seq_along(train_RT)) {
    
    # Load training data
    dt_RT <- train_RT[i] %>%
      vroom(col_names = c("peptide", "RT"), show_col_types = F, skip = 1)
    
    # Calibrate
    RCs = pyteomics$achrom$get_RCs(sequences = r_to_py(dt_RT$peptide), 
                                   RTs=r_to_py(dt_RT$RT), 
                                   term_aa=r_to_py(FALSE))
    # Save model parameters
    saveRDS(RCs, paste0(out_RT[i], "/achrom_RCs.rds"))
    rm(RCs)
  }
  
  # Predict
  for (i in seq_along(pred_RT)) {
    
    # Load parameters
    RCs <- readRDS(paste0(str_replace_all(pred_RT[i], 
                                          pattern = "results/RT_prediction/predict/", 
                                          replacement = "results/RT_prediction/RT_models/"), 
                          "/achrom_RCs.rds"))
    
    # Predict
    dt_test <- paste0(str_replace_all(pred_RT[i], pattern = "results/RT_prediction/predict/", replacement = "results/RT_prediction/test/"), ".tsv") %>%
      vroom(col_names = c("peptide", "RT"), show_col_types = F, skip = 1) 
    
    x <- dt_test %>%
      pull(peptide)  %>%
      r_to_py()
    
    py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
")
    dt_test$RT_pred = py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE)) 
    vroom_write(dt_test, file = paste0(pred_RT[i], "/test.tsv"), delim = "\t", num_threads = Ncpu)
    
  }
  # Empty cmds for Snakemake output rule
  cmd_RT_train <- NA
  cmd_RT_test <- NA
} 

### ---------------------------- (4) Export cmds --------------------------------------
dataset_name <- str_split_fixed(train_RT, pattern = "results/RT_prediction/", 2)[,2] %>%
  str_split_fixed(pattern = "/", 2)
dataset_name <- dataset_name[,2] %>%
  str_remove_all("/") %>%
  str_remove_all(".tsv")

if (method == "achrom") {
  dataset_name <- NA
}

# {
# data.frame(RT_dataset = dataset_name,
#            cmd = cmd_RT_train,
#            number = seq_along(cmd_RT_train))  %>%
#   vroom_write(delim = ",", append = FALSE,
#               file = "results/RT_prediction/cmd_RT_train.csv")
# 
# data.frame(number = seq_along(cmd_RT_test),
#            cmd = cmd_RT_test) %>%
#   vroom_write(delim = ",", append = FALSE,
#               file = "results/RT_prediction/cmd_RT_test.csv")
# }

data.frame(RT_dataset = dataset_name,
           cmd = cmd_RT_train,
           number = seq_along(cmd_RT_train)) %>%
  vroom_write(delim = ",", append = FALSE,
              file = unlist(snakemake@output[["cmd_RT_train"]]), col_names = T)

data.frame(number = seq_along(cmd_RT_test),
           cmd = cmd_RT_test)  %>%
  vroom_write(delim = ",", append = FALSE,
              file =  unlist(snakemake@output[["cmd_RT_test"]]), col_names = T)
