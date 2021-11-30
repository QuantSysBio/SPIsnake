### ---------------------------------------------- Define peptide aggregation  ----------------------------------------------
# description:  Evaluate RT prediction error. 
#               
# input:        1. RT command
#               2. RT Calibration peptides split into train/test
# output:       
#               A table with a single line per combination of parameters for peptide across proteome chunks. 
#               They will be used as wildcards to control uniqueness, MW computation, m/z matching and PTM generation. 
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
#suppressPackageStartupMessages(library(parallel))
require("bettermc")
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(vroom))
print(sessionInfo())

{
  ### Manual startup
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # cmd_RT_test <- vroom("results/RT_prediction/cmd_RT_test.csv", delim=',', show_col_types = FALSE)
  # method = "achrom"
}

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

# Cmd table
cmd_RT_test <- vroom(snakemake@input[["cmd_RT_test"]], delim=',', show_col_types = FALSE)

# Params
method = as.character(snakemake@params[["method"]])

RT_predictors <- c()
gg_out <- list()

### ---------------------------- (1) RT: Predict cmds ----------------------------
if (!is.na(cmd_RT_test$cmd)) {
  for (i in 1:length(cmd_RT_test$cmd)) {
    system(cmd_RT_test$cmd[i])
  }
}

### ---------------------------- (2) Performance evaluation --------------------------------------
# Ground truth
files_train <- list.files("results/RT_prediction/train/", pattern = ".tsv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(files_train) <- list.files("results/RT_prediction/train/", pattern = ".tsv", full.names = F) %>%
  str_remove(".tsv")

files_test <- list.files("results/RT_prediction/test/", pattern = ".tsv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(files_test) <- list.files("results/RT_prediction/test/", pattern = ".tsv", full.names = F) %>%
  str_remove(".tsv")

# Evaluate each sampling iteration separately
split_structure <- files_test %>%
  names() %>%
  as_tibble() %>%
  separate(value, sep = ".sample", into = c("dataset", "sample")) %>%
  mutate(sample = as.integer(str_remove(sample, "sample"))) %>%
  select(dataset, sample) %>%
  unique()

# RT
pred_RT <- list.files("results/RT_prediction/predict/", pattern = "test.tsv", full.names = T, recursive = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_RT) <- list.files("results/RT_prediction/predict/", pattern = "test.tsv", full.names = F, recursive = T) %>%
  str_remove_all(pattern = "/test.tsv")

if (unique(c("x", "y", "y_pred") %in% colnames(pred_RT[1]))) {
  # AutoRT input
  RT <- pred_RT %>%
    bind_rows(.id = "file") %>%
    dplyr::rename(peptide = x)  %>%
    mutate(RT = y/60,
           RT_pred = y_pred) %>%
    select(-RT) %>%
    split(~file)
  
} else if (unique(c("peptide", "RT", "RT_pred") %in% colnames(pred_RT[[1]]))) {
  RT <- pred_RT %>%
    bind_rows(.id = "file") %>%
    select(-RT) %>%
    split(~file)
}

RT_predictors <- c(RT_predictors, "RT")

out <- tibble()
for (i in seq_along(split_structure$dataset)) {
  
  dataset <- split_structure$dataset[i]
  keep <- files_test[grep(dataset, names(files_test))]
  
  for (j in seq_along(keep)) {
    
    for (pred in RT_predictors) {
      
      dt = keep[[j]] %>%
        left_join(get(pred)[[names(keep)[j]]], by="peptide") %>%
        na.omit()
      
      out_tmp <- regression_stats(obs = dt$RT, 
                                  pred = dt$RT_pred) %>%
        t() %>%
        as_tibble() %>%
        mutate(dataset = dataset,
               sample = names(keep)[j],
               predictor = method)
      out <- rbind(out, out_tmp)
    }
  }
}
out <- out %>%
  mutate(predictor = str_remove(predictor, "pred_"),
         sample = as.numeric(str_split_fixed(sample, ".sample", 2)[,2])) %>%
  pivot_longer(cols = c("Rsquared", "PCC", "MSE", "RMSE", "MAE"), names_to = "metric")

# Print
RT_Performance_df <- out %>% 
  unique() %>% 
  group_by(dataset, metric) %>% 
  summarise(mean_value = mean(value))
RT_Performance_df %>% print.data.frame()

metrics = c("Rsquared", "PCC", "MSE", "RMSE", "MAE")
gg <- list()
for (i in seq_along(metrics)) {
  
  keep_col = metrics[[i]]
  gg[[i]] <- out %>%
    filter(metric == keep_col) %>%
    ggplot(aes(y=value, x=as.factor(dataset))) + 
    geom_boxplot(aes(fill=as.factor(predictor)), stat="boxplot", position="dodge", alpha=0.5, width=0.2) + 
    theme_bw() + 
    theme(text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
    guides(fill=guide_legend(title="predictor")) + 
    xlab("dataset") + 
    ylab(keep_col)
  names(gg)[i] = keep_col
  
}
gg_out$performance <- arrangeGrob(grobs = gg[c(1,2,5)], ncol = 1)

### Scatterplots
RT_predictors <- c("files_test", RT_predictors)

dt <- lapply(RT_predictors, get) %>%
  lapply(bind_rows, .id = "file") 
names(dt) <- RT_predictors
dt <- dt[2:length(dt)] %>%
  bind_rows(.id = "method") %>%
  left_join(dt$files_test) %>%
  mutate(dataset = str_split_fixed(file, pattern = fixed("."), n = 3)[,1]) %>%
  mutate(sample = str_split_fixed(file, pattern = fixed("."), n = 3)[,2]) %>%
  select(-file)

gg_out$Scatter_pred_obs <- ggplot(dt, aes(y=RT_pred, x=RT)) + 
  geom_point(stat="identity", position="identity", alpha=0.3, size=1, color="firebrick") + 
  geom_smooth(stat="smooth", position="identity", method="lm", se=TRUE, n=80, level=0.95, span=0.75, alpha=1) + 
  facet_grid(dataset ~ method) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", size=11, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 1)) + 
  ggtitle("Retention Time: Predicted vs Observed") + 
  xlab("RT") + 
  ylab("RT_pred")

### ---------------------------- (3) Export --------------------------------------
# Performance plots
for (i in 1:length(gg_out)) {
  print(names(gg_out)[i])
  ggsave(plot = gg_out[[i]], filename = paste0("results/RT_prediction/plots/",names(gg_out)[i],".pdf"), width = 16, height = 9, dpi = "retina")
  ggsave(plot = gg_out[[i]], filename = paste0("results/RT_prediction/plots/",names(gg_out)[i],".png"), width = 16, height = 9, dpi = "retina")
}

# Performance dt
vroom_write(RT_Performance_df, delim = ",", append = FALSE,
            file = unlist(snakemake@output[["RT_Performance_df"]]))
