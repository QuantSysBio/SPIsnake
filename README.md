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
```
## To do:
- [x] Proteome preprocessing: split sequences by max length
- [x] Proteome sorting with Linclust
- [x] Chunk generation
- [ ] Control peptide generation rule via Master_sheet (proteome and corresponding peptide lengths) and checkpoint over chunks
- [ ] Peptide generator
- - [x] Paralell PCP/PSP generation
- - [X] PCP/PSP uniqueness inside chunks and save .csv.gz using AA index
- - [ ] PCP/PSP control in peptide generator
- - [ ] Gather generated peptides across chunks and save unique peptide sequences
- [ ] MW computation
- [ ] PTM sequence generation (control via Master_sheet)
- [ ] MW filtering
- [ ] RT prediction
- [ ] IC50 prediction
