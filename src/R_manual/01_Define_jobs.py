# -*- coding: utf-8 -*-
"""
Temp script to work with objects in Python

"""

import glob
import pandas as pd
import os, multiprocessing
import yaml
import snakemake
from snakemake.utils import min_version


#%%
config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
runfile(snakefiles + "01_folders.py")
#include: snakefiles + "02_DB_wildcards.py"
#include: snakefiles + "03_generate_exhaustive_DB.py"
#include: snakefiles + "04_filter_DB.py"


#%%  

def Expand_df_column(filepath, col_name, prefix, extension) :
    df = pd.read_csv(filepath)
    
    row = []
    for i in range(len(df)):
        suffix = df.loc[i, col_name]
        row += [prefix + suffix + extension]
        
    return(row)

tmp = Expand_df_column(filepath = "/home/yhorokh/Snakemake/SPI-snake/Master_table.csv", 
                   col_name = "Proteome",
                   prefix = "/home/yhorokh/Snakemake/SPI-snake/data/reference/",
                   extension = ".fasta")


#%%  

Master_table = (
    pd.read_csv(features["Master_table"], sep=",")
)

def get_proteomes_input(Master_table) :
    proteomes = expand("data/reference/{Proteome}.fasta", Proteome=Master_table["Proteome"])
    return(proteomes)

print(get_proteomes_input(Master_table))




#%%
def get_Master_table():
    Master_table_expanded = (
    pd.read_csv(os.path.join(dir_DB_exhaustive, "Master_table_expanded.csv"), 
        sep=",", 
        dtype={
        "Proteome": str, 
        "Splice_type": str, 
        "N_mers": int, 
        "Min_Interv_length": int, 
        "chunk": str})
    .set_index(["Proteome", "Splice_type", "N_mers", "Min_Interv_length", "chunk"], drop=False)
    .sort_index()
    )
    return Master_table_expanded

MasterTable = get_Master_table()
#MT = MasterTable.loc[Proteome].loc[Splice_type].loc[N_mers].loc[Min_Interv_length].loc[chunk]

def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(MasterTable.loc[wildcards.Proteome]))
    )
    assert len(files) > 0
    return files


Master_table_expanded = (
pd.read_csv(os.path.join(dir_DB_exhaustive, "Master_table_expanded.csv"), 
    sep=",", 
    dtype={
    "chunk": str}).set_index(["chunk"], drop=False).sort_index()
)

Master_table_expanded = pd.read_csv(os.path.join(dir_DB_exhaustive, "Master_table_expanded.csv"), 
    sep=",", 
    dtype={
    "Proteome": str, 
    "Splice_type": str, 
    "N_mers": int, 
    "Min_Interv_length": int, 
    "chunk": str}).set_index(["Proteome", "Splice_type", "N_mers", "Min_Interv_length", "chunk"], drop=False).sort_index()

#%%


pattern = expand(Proteome=Proteome, Splice_type=Splice_type, N_mers=N_mers, Min_Interv_length=Min_Interv_length, chunk=chunk, **w).format(
                Proteome=MT.Proteome, Splice_type=MT.Splice_type, N_mers=MT.N_mers, Min_Interv_length=MT.Min_Interv_length, chunk=MT.chunk,
            )
tmp = join(dir_DB_exhaustive, pattern)

print(tmp)


        
os.path.join()
