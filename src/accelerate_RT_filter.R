library(profvis)

# MW computation
profvis({
  PSP <- PSP_list %>%
    rbindlist() %>%
    lazy_dt() %>%
    unique() %>%
    mutate(MW=computeMZ_biostrings(peptide)) %>%
    as.data.table()
})

profvis({
  PSP <- PSP_list %>%
    rbindlist() %>%
    mutate(index = str_sub(peptide, index_length + 1, index_length + 2)) %>%
    as.data.table() %>%
    split(by = c("index"), drop = T) %>%
    mclapply(mc.cores = Ncpu, FUN = function(x){
      x %>%
        lazy_dt() %>%
        select(-index) %>%
        unique() %>%
        mutate(MW=computeMZ_biostrings(peptide)) %>%
        as.data.table()
    }) %>%
    rbindlist()
})

# ----------------------------------------------- MW filter --------------------------------------------
profvis({
  x <- input[MW %inrange% mzList[,c("MW_Min", "MW_Max")]]
})

profvis({
  x <- input %>%
    lazy_dt() %>%
    mutate(index = str_sub(peptide, index_length + 1, index_length + 2)) %>%
    as.data.table() %>%
    split(by = c("index"), drop = T) %>%
    mclapply(mc.cores = Ncpu, FUN = function(x){
      x = x %>%
        lazy_dt() %>%
        select(-index) %>%
        as.data.table()
      x = x[MW %inrange% mzList[,c("MW_Min", "MW_Max")]]
      return(x)
    }) %>%
    rbindlist()
})


# ----------------------------------------------- RT compute --------------------------------------------
# Predict
py_calls <- py_run_string("
def achrom_calculate_RT(x, RCs, raise_no_mod):
  x = pd.DataFrame({'sequences': x})
  out = x['sequences'].apply(
    lambda x : achrom.calculate_RT(x, RCs, raise_no_mod=False)
  )
  return out
")
profvis({
y <- mz_nomod[[i]][[MS_mass_list]] %>%
  lazy_dt() %>%
  mutate(index = str_sub(peptide, index_length + 1, index_length + 2)) %>%
  as.data.table() %>%
  split(by = c("index"), drop = T) %>%
  mclapply(mc.cores = Ncpu, FUN = function(x){
    x = x %>%
      lazy_dt() %>%
      select(-index) %>%
      pull(peptide) %>%
      r_to_py()
    x = data.table(RT = py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE)))
    return(x)
  }) %>%
  rbindlist() %>%
  pull(RT)
})

profvis({
  x <- mz_nomod[[i]][[MS_mass_list]] %>%
    pull(peptide) %>%
    r_to_py()
  y = py_calls$achrom_calculate_RT(x, RCs, r_to_py(FALSE))
})




# ----------------------------------------------- 2D filter --------------------------------------------
profvis({
  y <- mz_nomod[[i]][[MS_mass_list]][
    MW %inrange% mzList[,c("MW_Min", "MW_Max")] & RT_pred %inrange% mzList[,c("RT_Min", "RT_Max")]]
})

profvis({
  y2 <- mz_nomod[[i]][[MS_mass_list]] %>%
    lazy_dt() %>%
    mutate(index = str_sub(peptide, index_length + 1, index_length + 2)) %>%
    as.data.table() %>%
    split(by = c("index"), drop = T) %>%
    mclapply(mc.cores = Ncpu, FUN = function(x){
      x = x %>%
        lazy_dt() %>%
        select(-index) %>%
        as.data.table()
      x = x[MW %inrange% mzList[,c("MW_Min", "MW_Max")] & RT_pred %inrange% mzList[,c("RT_Min", "RT_Max")]]
      return(x)
    }) %>%
    rbindlist()
})
