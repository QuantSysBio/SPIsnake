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



Master_table = (
    pd.read_csv(features["Master_table"], sep=",", dtype={"Proteome": str})
    .set_index("Proteome", drop=False)
    .sort_index()
)



#%%  

def get_proteomes_input(wildcards) :
    Proteome = Master_table.loc[wildcards.Proteome]
    return expand("data/reference/{Proteome}.fasta", Proteome=Proteome)





#%%
tmp = Expand_df_column(filepath = "/home/yhorokh/Snakemake/SPI-snake/Master_table.csv", 
                   col_name = "Proteome",
                   prefix = "/home/yhorokh/Snakemake/SPI-snake/data/reference/",
                   extension = ".fasta")


