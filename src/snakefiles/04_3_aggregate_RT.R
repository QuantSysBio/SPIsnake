### ---------------------------------------------- Define peptide aggregation  ----------------------------------------------
# description:  Evaluate RT prediction error. 
#               
# input:        1. AutoRT command
#               2. RT Calibration peptides split into train/test
# output:       
#               A table with a single line per combination of parameters for peptide across proteome chunks. 
#               They will be used as wildcards to control uniqueness, MW computation, m/z matching and PTM generation. 
#               
# author:       YH

### Log
log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(parallelly))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(data.table))

{
  ### Manual startup
  # setwd("/home/yhorokh/SNAKEMAKE/SPIsnake")
  # cmd_AutoRT_test <- vroom("results/RT_prediction/cmd_AutoRT_test.csv", delim=',')
}

source("src/snakefiles/functions.R")
print("Loaded functions. Loading the data")

# Cmd table
cmd_AutoRT_test <- vroom(snakemake@input[["cmd_AutoRT_test"]], delim=',')

### ---------------------------- (1) AutoRT: Predict cmds ----------------------------
for (i in 1:length(cmd_AutoRT_test$cmd)) {
  system(cmd_AutoRT_test$cmd[i])
}

### ---------------------------- (2) Performance evaluation --------------------------------------
# Ground truth
files_test <- list.files("results_MeV/test/achrom/", pattern = ".csv", full.names = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(files_test) <- list.files("results_MeV/test/achrom/", pattern = ".csv", full.names = F) %>%
  str_remove(".csv")



# AutoRT
pred_AutoRT <- list.files("results_MeV/predict/AutoRT/", pattern = "test.tsv", full.names = T, recursive = T) %>%
  lapply(vroom, show_col_types = FALSE)
names(pred_AutoRT) <- list.files("results_MeV/predict/AutoRT/", pattern = "test.tsv", full.names = F, recursive = T) %>%
  str_remove_all(pattern = "/test.tsv")

pred <- lapply(pred_AutoRT, function(x) predict(lm(y ~ y_pred, data = x))) %>%
  lapply(as_tibble) %>%
  bind_rows()

pred_AutoRT <- pred_AutoRT %>%
  bind_rows(.id = "file") %>%
  rename(Mascot_seq = x)  %>%
  mutate(RT_pred = pred$value) %>%
  select(-y) %>%
  select(-y_pred) %>%
  split(~file)

RT_predictors <- c(RT_predictors, "AutoRT")


out <- tibble()
for (i in seq_along(split_structure$dataset)) {
  
  dataset <- split_structure$dataset[i]
  keep <- files_test[grep(dataset, names(files_test))]
  
  for (j in seq_along(keep)) {
    
    for (pred in RT_predictors) {
      
      dt = keep[[j]] %>%
        left_join(get(pred)[[names(keep)[j]]]) %>%
        na.omit()
      
      out_tmp <- regression_stats(obs = dt$RT, 
                                  pred = dt$RT_pred) %>%
        t() %>%
        as_tibble() %>%
        mutate(dataset = dataset,
               sample = names(keep)[j],
               predictor = pred)
      out <- rbind(out, out_tmp)
    }
  }
}
out <- out %>%
  mutate(predictor = str_remove(predictor, "pred_"),
         sample = str_split_fixed(sample, ".sample", 2)[,2]) %>%
  pivot_longer(cols = c("Rsquared", "PCC", "MSE", "RMSE", "MAE"), names_to = "metric")


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
gg[1]
gg[2]
gg[3]
gg[4]
gg[5]

gg$performance <- arrangeGrob(grobs = gg[c(1,2,5)], ncol = 1)
plot(gg$performance)

### Scatterplots
RT_predictors <- c("files_test", RT_predictors)
RT_predictors

dt <- lapply(RT_predictors, get) %>%
  lapply(bind_rows, .id = "file") 
names(dt) <- RT_predictors
dt <- dt[2:length(dt)] %>%
  bind_rows(.id = "method") %>%
  left_join(dt$files_test) %>%
  mutate(dataset = str_split_fixed(file, pattern = fixed("."), n = 3)[,1]) %>%
  mutate(sample = str_split_fixed(file, pattern = fixed("."), n = 3)[,2]) %>%
  select(-file)
dt

library(ggpubr)
gg$scatter_correlation <- ggscatter(dt, x="RT", y="RT_pred", 
                  fill = "lightgrey", 
                  color = "black", shape = 21, size = 2, alpha = 0.3,
                  add.params = list(color = "orange", fill = "blue"), # Customize reg. line
                  facet.by = c("dataset", "method"),
                  title = paste("Predicted vs Observed retention time"), 
                  conf.int = T, 
                  conf.int.level = 1 - 10^-15,
                  # conf.int.level = 0.95,
                  add = "reg.line",
                  cor.coef = T, 
                  ggtheme = theme_bw())
plot(gg$scatter_correlation)

for (i in 1:5) {
  ggsave(gg[[i]], filename = paste0("results_MeV/plots/",names(gg)[i],".pdf"), width = 8, height = 6, dpi = "retina")
  ggsave(gg[[i]], filename = paste0("results_MeV/plots/",names(gg)[i],".png"), width = 8, height = 6, dpi = "retina")
}

for (i in 6:7) {
  ggsave(gg[[i]], filename = paste0("results_MeV/plots/",names(gg)[i],".pdf"), width = 16, height = 9, dpi = "retina")
  ggsave(gg[[i]], filename = paste0("results_MeV/plots/",names(gg)[i],".png"), width = 16, height = 9, dpi = "retina")
}


### ---------------------------- (3) Export --------------------------------------
# Performance dt
vroom_write(Peptide_aggregation_table, delim = ",", append = FALSE,
            file = unlist(snakemake@output[["RT_Performance_df"]]))

# Performance plots


# rule aggregate_RT:
#   input:
#   checkpoint = Checkpoint_RT_train(join(dir_RT_prediction, "AutoRT_models/{RT_dataset}/best_model.hdf5")),
# cmd_AutoRT_test = join(dir_RT_prediction, "cmd_AutoRT_test.csv")
# output:
#   RT_Performance_plot = join(dir_RT_prediction, "RT_Performance.pdf"),
# RT_Performance_df = join(dir_RT_prediction, "RT_Performance.csv")
# benchmark: 
#   join(benchmarks, "aggregate_RT.json")
# log: 
#   join(logs, "aggregate_RT.txt")
# conda: 
#   "R_env_reticulate.yaml"
# resources: # 1 per node at the time
#   load = 100 
# script:
#   "04_3_aggregate_RT.R"

# vroom_write(Peptide_aggregation_table, delim = ",", append = FALSE,
#             file = paste0(dir_DB_PTM_mz, "/Peptide_aggregation_table.csv"))
