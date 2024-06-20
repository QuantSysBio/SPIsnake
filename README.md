# SPIsnake
#### Spliced Peptide Identification, Search space Navigation And K-mer filtering Engine
SPIsnake performs an in-silico generation of the peptide search space and filters it in a data-driven manner to provide users with FASTA files, peptide annotation and statistics of the search space. It is intended for the application upstream of the MS search engines to explore and pre-filter the MS search space.

## Installation
```
mamba create -c conda-forge -c bioconda -n SPIsnake snakemake=7.25.4 python=3.10.10
mamba activate SPIsnake
mamba install -c conda-forge singularity=3.8.6
```

## Setup
1. Define sequence search space in `Master_table.csv`.  
N-mers separated by underscores will give rise to a full range of peptide sequences. Redundancies are removed internally.  
Place proteome `.fasta` files in `data/reference`

2. Specify the MS datasets and corresponding tolerances to be used for MW filtering in `data/Experiment_design.csv`; 
If the RAW files are available, place the in `data/RAW` and provide the PEAKS database search output to `data/DB search psm.csv`.  
Alternatively, place the `.txt` files with observed MS1 masses and peptide-RT pairs into `data/MS_mass_lists` and `data/RT_calibration`.  
If the empty placeholder file `data/DB search psm.csv` is replaced, the table contents will be exported to `data/RT_calibration`.

4. Define common job parameters in `features.yaml`
5. Specify the resource usage in `config.yml`
6. (optional) Place PTM tables in `.csv` format to `data/modifications` and  `data/modifications_fixed`. The names must match the entries in `data/Experiment_design.csv` and `Master_table.csv`.
7. (optional) Slurm parameters: update the `src/cluster.yaml`. Make sure that `output` directory exists and that `chdir` is the workflow directory.

## Execution
### Workflow check
`snakemake --use-singularity -j 1 -n`

### Workflow execution (local):
```
time snakemake --use-singularity -j 1 --resources load=100
```
### Workflow execution (Slurm):
```
time snakemake --use-singularity --cluster-config src/cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} --mem {cluster.mem} --job-name {cluster.job-name} -o {cluster.output} -D {cluster.chdir}" -j 3 -w 60 --restart-times 1 --cluster-cancel "scancel" --resources load=100
```

## Detailed input description
### Master table
| Parameter    | Meaning    | 
| ----------- | ----------- | 
| Proteome | FASTA file name in `data/reference` | 
| enzyme_type | The rule to digest the proteins. Use `PCP` for nonspecific cleavage, `PSP` or `cis-PSP` for cis-splicing. Other enzymes are called according to "cleaver" R package | 
| N_mers | Either a single integer or a pair of underscore-separated integers specifying the range of peptide lengths | 
| Max_Interv_length | Maximal length between splice reactant termini during the cis-splicing. In addition, used as a grouping variable in the final arrow dataset. This value must be non-empty. |
| MaxE | Hyperparameter controling the size of internal chunks for dataset processing |
| PTMs | Variable PTM set to be used for a given FASTA file |

### Experiment design
| Parameter    | Meaning    | 
| ----------- | ----------- | 
| Filename | RAW file name | 
| Biological_group | (Optional) A grouping variable to create a common output across multiple files | 
| Precursor_mass_tolerance_ppm | Numeric | 
| MHC-I_alleles | NetMHCPan allele names. Multiple alleles can be provided with pipe separator | 
| Affinity_threshold | A threshold for predicted MHC-I affinity (nM). Multiple alleles can be provided with pipe separator | 
| modifications_fixed | Fixed PTM set to be used for a given data file | 

### features
| Parameter    | Meaning    | 
| ----------- | ----------- | 
| min_protein_length | Minimal length of protein to process | 
| max_protein_length | Maximal length of protein to process. Longer proteins will be internally split into fixed length chunks | 
| replace_I_with_L | Whether to treat Ile as Leu (not recommended for large searches) | 
| max_variable_PTM | Maximal number of variable PTMs per peptide | 
| generate_spliced_PTMs | Whether to generate variable PTMs on spliced peptides (not recommended for large searches) | 
| PTM_chunk | Batch size when generating PTMs | 
| netMHCpan_path | Location of NetMHCpan inside the Docker container | 
| netMHCpan_chunk | Batch size when predicting MHC-I affinity | 
| AA_index_length | The number of amino acids to use for database partitioning. It's recommended to use 1 for enzymatic cleavage and 2 for large inputs or cis-splicing | 
| duckdb_RAM | Fraction of available RAM to use by duckDB | 
| duckdb_max_retries | Number of reties for executing duckDB aggregations | 
| ThermoRawFileParser | Location of ThermoRawFileParser inside the Docker container | 
| method | Whether to use `achrom` or `AutoRT` for retention time prediction | 
| n_folds | Number of folds for RT predictor cross-validation | 
| proportion_train | Fraction of peptide sequences in the train set for RT predictor | 
| quantile | Stringency of RT filter cutoff. The quantile of RT prediction error | 
| FASTA_outputs | Whether to generate a FASTA file for a given filtering step | 
| strata_sizes | Whether to estimate the size of unfiltered search spaces | 
| strata_sizes | Whether to estimate the size of filtered search spaces | 
| filtering_sizes | Whether to estimate the size of unfiltered search spaces | 

### config
| Parameter    | Meaning    | 
| ----------- | ----------- | 
| max_cpus | Maximal number of available CPUs to use | 
| cpus_critical | Keep at 1 | 
| max_mem | Maximal amount of RAM to use  | 
| min_mem | Minimal amount of RAM to use | 

## Outputs
The peptide-level filtering information is stored as a Hive-partitioned arrow dataset at: `results/DB_out/arrow/`
Every peptide is associated with molecular weight (MW) and retention time (RT), (optional) MHC-I affinity and whether the peptide sequence could be explained by the MS data. 

The protein-peptide mapping is stored as a Hive-partitioned arrow dataset at: `results/DB_out/peptide_mapping/`

FASTA files: `results/DB_out/FASTA/`

Database aggragation stats: `results/DB_out/Stats/`
Reported is the number of unique peptides across the FASTA files and filtering steps.
