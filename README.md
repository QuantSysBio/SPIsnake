# SPIsnake


## Installation
```
conda create -c conda-forge -c bioconda -n SPIsnake snakemake
conda activate SPIsnake
conda install -c conda-forge singularity 
```
For cluster execution, please install mamba package manager: `conda install -n base -c conda-forge mamba`.

## Setup
1. Define common job parameters in `config.yml` and `features.yml`
2. Define proteomic search space in `Master_table.csv`. 
N-mers separated by underscores will give rise to a full range of peptide sequences. Potentially redundancies are removed internally. 
Place proteome `.fasta` files in `data/reference`
Place PTM tables `.csv` files in `data/modifications`
3. Specify datasets and corresponding tolerances to be used for m/z matching in `/data/Experiment_design.csv` and place the corresponding `.txt` files in `/data/MS_mass_lists`


## Execution
```
snakemake --use-conda -j 1 -n
time snakemake --use-conda -j 27 --conda-frontend conda --resources load=100
time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100
```

## SPIsnake + Slurm
Connect to the Mascot server first: `ssh username@s1604-mascot01`.
Slurm jobs are currently submitted from the fileserver. Therefore, enter: `ssh username@s1604-fs01` and navigate to `/data/SPIsnake`.

### Install conda in your home directory
Enter `bash Miniconda3-latest-Linux-x86_64.sh` and follow the instructions.
After that, create the SPIsnake environment as described under **Installation**.

### Clone repo + upload data
Enter `git clone https://github.com/QuantSysBio/SPIsnake` to retrieve the latest code.
If necessary, deposit your data using `sftp`.

### Jobscript
The `jobscript.sh` script specifies parameters of both Slurm and Snakemake.
More info about Slurm's `sbatch` can be found [here](https://slurm.schedmd.com/sbatch.html). Useful guidelines about creating jobscripts can be found [here](https://docs.gwdg.de/doku.php?id=en:services:application_services:high_performance_computing:courses:scc-introductory-course).
You can get an overview about which compute nodes are assigned to which partition by calling `sinfo`.

### Submit jobs
Make sure you have activated the SPIsnake environment.
Enter `sbatch jobscript.sh`. You can track the status of your job via the outfile of by `ssh`-ing on the compute nodes and monitoring their activity.
