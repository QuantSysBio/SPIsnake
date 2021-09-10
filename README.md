# SPIsnake


## Installation
```
conda create -c conda-forge -c bioconda -n SPIsnake snakemake
conda activate SPIsnake
conda install -c conda-forge singularity 
```

## Execution
```
snakemake --use-conda -j 1 -n
time snakemake --use-conda -j 27 --conda-frontend conda --resources load=100
time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100
```
