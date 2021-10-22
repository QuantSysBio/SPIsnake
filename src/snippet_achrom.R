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
