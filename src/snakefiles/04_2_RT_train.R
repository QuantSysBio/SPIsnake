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
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # RT_dataset = "results/RT_prediction/AutoRT_models/MeV_MA0009-BE08_allFractions/model.json"
  # cmd_AutoRT_train = vroom(file = "results/RT_prediction/cmd_AutoRT_train.csv", delim=',', show_col_types = FALSE)
}
source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

### ---------------------------- (1) Read inputs ----------------------------
# Wildcard
RT_dataset = snakemake@output[[1]]

# Cmd table
cmd_AutoRT_train <- vroom(snakemake@input[["cmd_AutoRT_train"]], delim=',', show_col_types = FALSE)

# Output dirs
dir_RT_calibration = snakemake@params[["dir_RT_calibration"]]

### ---------------------------- (2) Define aggregation wildcards --------------------------------------
# Choose the train cmd
RT_dataset <- str_split_fixed(RT_dataset, pattern = "results/RT_prediction/AutoRT_models/", n = 2)[,2] %>%
  str_remove_all(pattern = "/model.json") 

# Call AutoRT
system(cmd_AutoRT_train$cmd[cmd_AutoRT_train$RT_dataset == RT_dataset])
