shell.executable("/bin/bash")

#container: "docker://snakemake/snakemake:v6.11.1"
container: "docker://yhorokh/spisnake:3.0.0"
import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

min_version("6.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
include: snakefiles + "01_folders.py"
include: snakefiles + "02_DB_wildcards.py"
include: snakefiles + "03_RT_prediction.py"

rule all:
    input:
        join(dir_DB_exhaustive, ".Generate_peptides.done"),
        join(dir_DB_out, ".Aggregare_FASTA.done"),
        join(dir_DB_out, ".Aggregare_peptide_mapping.done"),
        join(dir_DB_out, ".Aggregare_Stats.done")


### Execution
# snakemake -j 1 -n
# time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100

### Plots
# snakemake --report
# snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
# snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg
# rm filegraph.dot
# rm dag.dot
