### ---------------------------------------------- Save outputs ----------------------------------------------
# description:  Save SPIsnake main outputs
# 
# output:       Arrow dataset; peptide sequences with filter info
#               
# author:       YH

### ---------------------------- (1) Peptde information: arrow --------------------------------------
### Peptide info
cat(as.character(Sys.time()), " - ", "Start saving peptide sequences with filter info", "\n")
if ("chunks" %in% colnames(pep)) {
  pep$chunks <- NULL
}

pep <- pep[, c("proteome", "MiSl", "chunk") := list(chunk_params$proteome, chunk_params$MiSl, chunk_params$chunk)] 
setcolorder(pep, neworder = c("index", "length", "proteome", "enzyme", "MiSl", "chunk", "peptide", "MW.unmodified"))
pep[1:nrow(pep), write_dataset(group_by(as_arrow_table(.SD), 
                                        length, proteome, enzyme, MiSl, chunk), 
                               path = paste0(dir_DB_PTM_mz,  "/peptide_seqences/index=", .BY), 
                               existing_data_behavior = "overwrite",
                               format = "parquet", 
                               max_partitions = 10240L,
                               compression = "lz4"), 
    by=index, 
    .SDcols=colnames(pep)]
cat(as.character(Sys.time()), " - ", "Saved peptide sequences with filter info: Done\n")

### Clean-up
suppressWarnings(rm(pep_bg, fa, PTMs, bg, Biological_group, Biological_group_bg, select_mass_list))
