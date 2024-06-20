### ---------------------------------------------- RT train  ----------------------------------------------
# description:  Launch AutoRT train
#               
# input:        1. AutoRT command
#               2. RT Calibration peptides split into train/test
# output:       
#               - AutoRT model to be used for prediction
#               
# author:       YH

### Log
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, append=TRUE, type=c("output", "message"), split = TRUE)
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Read inputs ----------------------------
if (exists("snakemake")) {
  # Wildcard
  RT_dataset = snakemake@output[[1]]
  out <- RT_dataset
  
  # Cmd table
  cmd_RT_train <- fread(snakemake@input[["cmd_RT_train"]], sep=',')
  
  # Output dirs
  dir_RT_calibration = snakemake@params[["dir_RT_calibration"]]
  
} else {
  ### Manual startup
  RT_dataset = "results/RT_prediction/RT_models/nan/model.json"
  out <- RT_dataset
  cmd_RT_train = fread(file = "results/RT_prediction/cmd_RT_train.csv", sep=',')
}

### ---------------------------- (2) Define aggregation wildcards --------------------------------------
# Choose the train cmd
RT_dataset <- str_split_fixed(RT_dataset, pattern = "results/RT_prediction/RT_models/", n = 2)[,2] %>%
  str_remove_all(pattern = "/model.json") 

if (!RT_dataset == "nan") {
  # Call external training cmd
  system(cmd_RT_train$cmd[cmd_RT_train$RT_dataset == RT_dataset])
  
} else {
  # Empty output for achrom filter
  suppressWarnings(dir.create(str_remove_all(out, pattern = "/model.json") ))
  system(paste("touch", out))
}
if (exists("snakemake")) {
  sink()
}
