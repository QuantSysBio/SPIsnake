# SPI-snake


## Installation
```
conda create -c conda-forge -c bioconda -n SPIsnake snakemake
conda activate SPIsnake
conda install -c conda-forge singularity 
```

## Execution
```
snakemake --use-conda -j 1
time snakemake make_all_files --use-conda -j 27 --conda-frontend conda --resources load=100
```
