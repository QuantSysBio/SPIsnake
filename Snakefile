shell.executable("/bin/bash")

singularity: "docker://continuumio/miniconda3:4.4.10"
import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

min_version("5.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
include: snakefiles + "00_functions.py"
include: snakefiles + "01_folders.py"
include: snakefiles + "02_DB_wildcards.py"
#include: snakefiles + "03_generate_exhaustive_DB.py"
#include: snakefiles + "04_filter_DB.py"

rule all:
    input:
        join(dir_DB_exhaustive, ".Split_proteomes.done")
        #Checkpoint_MakePattern(join(dir_DB_exhaustive, "{proteome}/{chunk}.fasta"))
        #Checkpoint_MakePattern(dir_DB_exhaustive + "inputSequence_{prot_chunk}.fasta")



### snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
### snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg


### snakemake --use-conda --use-singularity -r --verbose
### snakemake --filegraph | dot | display
### snakemake --dag | dot | display
