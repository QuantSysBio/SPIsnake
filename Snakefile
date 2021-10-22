shell.executable("/bin/bash")

container: "docker://snakemake/snakemake:v6.8.0"
import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

min_version("6.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
include: snakefiles + "00_functions.py"
include: snakefiles + "01_folders.py"
include: snakefiles + "02_DB_wildcards.py"
include: snakefiles + "03_PTM_mz_matching.py"
include: snakefiles + "04_RT_prediction.py"

rule all:
    input:
        join(dir_DB_exhaustive, ".Generate_indices.done"),
        join(dir_DB_PTM_mz, ".Aggregate_peptides.done"),
        join(dir_DB_exhaustive, ".Generate_peptides.done")


### Execution
# snakemake -j 1 -n
# time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100

### Plots
# snakemake --report
# snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
# snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg
# rm filegraph.dot
# rm dag.dot
