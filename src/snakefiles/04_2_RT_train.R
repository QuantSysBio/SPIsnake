### ---------------------------------------------- RT train  ----------------------------------------------
# description:  Launch AutoRT train
#               
# input:        1. AutoRT command
#               2. RT Calibration peptides split into train/test
# output:       
#               AutoRT model to be used for prediction
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log, append=TRUE, type=c("output", "message"))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vroom))

{
  #setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  #RT_dataset = "results/RT_prediction/RT_models/nan/model.json"
  #cmd_RT_train = vroom(file = "results/RT_prediction/cmd_RT_train.csv", delim=',', show_col_types = FALSE)
}
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")
print(sessionInfo())

### ---------------------------- (1) Read inputs ----------------------------
# Wildcard
RT_dataset = snakemake@output[[1]]
out <- RT_dataset

# Cmd table
cmd_RT_train <- vroom(snakemake@input[["cmd_RT_train"]], delim=',', show_col_types = FALSE)

# Output dirs
dir_RT_calibration = snakemake@params[["dir_RT_calibration"]]

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
