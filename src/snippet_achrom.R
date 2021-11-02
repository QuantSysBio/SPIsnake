dt_RT <- paste0("data/RT_calibration/", RT_calibration_lists$RT_list_file[1]) %>%
  vroom(col_names = c("row", "peptide", "RT", "scans"), show_col_types = F, skip = 1)

library(reticulate)
use_condaenv("r-reticulate")
pyteomics <- import("pyteomics")

py_run_string("
import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None
")

# Calibrate
RCs = pyteomics$achrom$get_RCs(sequences = r_to_py(dt_RT$peptide), 
                               RTs=r_to_py(dt_RT$RT), 
                               term_aa=r_to_py(FALSE))

# Predict
x = mz_nomod[[i]][[MS_mass_list]]$peptide %>%
  str_replace_all(pattern = "I", replacement = "L") %>%
  r_to_py()

py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
")
RT_pred = py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE)) %>% as.vector()


# Plot RCs vs achrom RCs
{
  df_RCs <- list()
  df_RCs$RCs_browne_hfba <- pyteomics$achrom$RCs_browne_hfba$aa 
  df_RCs$RCs_browne_tfa <- pyteomics$achrom$RCs_browne_tfa$aa 
  df_RCs$RCs_gilar_atlantis_ph10_0 <- pyteomics$achrom$RCs_gilar_atlantis_ph10_0$aa 
  df_RCs$RCs_gilar_atlantis_ph3_0 <- pyteomics$achrom$RCs_gilar_atlantis_ph3_0$aa 
  df_RCs$RCs_gilar_atlantis_ph4_5 <- pyteomics$achrom$RCs_gilar_atlantis_ph4_5$aa 
  df_RCs$RCs_gilar_beh <- pyteomics$achrom$RCs_gilar_beh$aa 
  df_RCs$RCs_gilar_beh_amide <- pyteomics$achrom$RCs_gilar_beh_amide$aa 
  df_RCs$RCs_gilar_rp <- pyteomics$achrom$RCs_gilar_rp$aa 
  df_RCs$RCs_guo_ph2_0 <- pyteomics$achrom$RCs_guo_ph2_0$aa 
  df_RCs$RCs_guo_ph7_0 <- pyteomics$achrom$RCs_guo_ph7_0$aa 
  df_RCs$RCs_krokhin_100A_fa <- pyteomics$achrom$RCs_krokhin_100A_fa$aa 
  df_RCs$RCs_krokhin_100A_tfa <- pyteomics$achrom$RCs_krokhin_100A_tfa$aa 
  df_RCs$RCs_meek_ph2_1 <- pyteomics$achrom$RCs_meek_ph2_1$aa 
  df_RCs$RCs_meek_ph7_4 <- pyteomics$achrom$RCs_meek_ph7_4$aa 
  df_RCs$RCs_palmblad <- pyteomics$achrom$RCs_palmblad$aa 
  df_RCs$RCs_yoshida <- pyteomics$achrom$RCs_yoshida$aa 
  df_RCs$RCs_yoshida_lc <- pyteomics$achrom$RCs_yoshida_lc$aa 
  df_RCs$RCs_zubarev <- pyteomics$achrom$RCs_zubarev$aa 
  # df_RCs$PCP_calibration <- RCs$aa

  dt_RCs <- df_RCs %>%
    lapply(function(x){
      x %>%
        as.data.table() 
    }) %>%
    bind_rows(.id = "achrom_dataset") %>%
    split(by=c("achrom_dataset")) %>%
    Reduce(function(...) merge(..., all = TRUE), x = .) %>% 
    as.data.frame()
  dist_RCs <- dt_RCs %>%
    select(-achrom_dataset)
  rownames(dist_RCs) <- dt_RCs$achrom_dataset
  dist_RCs <- dist_RCs %>%
    as.matrix() 
  
  pheatmap(dist_RCs, 
           na_col = "grey", cluster_rows = hclust(dist(dist_RCs)),
           labels_row=dt_RCs$achrom_dataset, 
           main = "Euclidean dist between RC datasets in AA parameter space")
}