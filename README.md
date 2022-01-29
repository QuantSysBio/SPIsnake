# SPIsnake


## Installation
```
conda create -c conda-forge -c bioconda -n SPIsnake snakemake
conda activate SPIsnake
conda install -c conda-forge singularity 
```

## Setup
1. Define common job parameters in `config.yml` and `features.yaml`
2. Define proteomic search space in `Master_table.csv`. 
N-mers separated by underscores will give rise to a full range of peptide sequences. Redundancies are removed internally. 
Place proteome `.fasta` files in `data/reference`
Place PTM tables `.csv` files in `data/modifications`
3. Specify datasets and corresponding tolerances to be used for MW filtering in `/data/Experiment_design.csv`; place the corresponding `.txt` files in `/data/MS_mass_lists`

## Execution
### Workflow check
`snakemake -j 1 -n`
### Workflow execution:
```
time snakemake --use-singularity -j 7 --resources load=100
```

## SPIsnake + Slurm
Connect to the Mascot server first: `ssh username@s1604-mascot01`.
Slurm jobs are currently submitted from the fileserver. Therefore, enter: `ssh username@s1604-fs01` and navigate to `/data/SPIsnake`.

### Install conda in your home directory
Enter `bash Miniconda3-latest-Linux-x86_64.sh` and follow the instructions.
After that, create the SPIsnake environment as described under **Installation**.

### Clone repo + upload data
Enter `git clone https://github.com/QuantSysBio/SPIsnake` to retrieve the latest code. You might need to [generate a token](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/) since GitHub recently removed the password authentication.
If necessary, deposit your data in the correct directory using `sftp`, `scp` or `rsync`. Instructions can be found in the [QSB getting started](https://pad.gwdg.de/s/JlkAOXJ2f#) document.

### Cluster execution
- Make sure you are in the correct directory (`data/SPIsnake/SPIsnake`) and on the correct node (`s1604-fs01`)
- Snakemake is executed from a Bash screen session that prevents the job from crashing once you disconnect from `ssh`. Therefore, enter:
`screen -S spisnake`
- Activate the conda environment:
`conda activate SPIsnake`
- Submit the job to the `elbe` partition. (You can get an overview about which compute nodes are assigned to which partition by calling `sinfo`.)  
```
snakemake --use-singularity --cluster-config src/cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} --mem {cluster.mem} --job-name {cluster.job-name} -o {cluster.output} -D {cluster.chdir} --exclusive" --conda-frontend conda -j 3 -w 36000
```
- Detach from the screen session by pressing `Ctrl+a+d`. You can resume to the session to check the progress via `screen -r spisnake`

### troubleshooting
In case you encounter an error with singularity that says `Write on output file failed because No space left on device`, go to the `/tmp` directory on `s1604-fs01`. Remove the directory `conda/` **after making sure that is was created by yourself and modified recently**. Then, try again.

Check the resource usage of current Slurm jobs via:
```
sacct --format='JobID,JobName,State,Elapsed,AllocNodes,NCPUS,NodeList,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask,ReqMem,MaxDiskWrite'
```

In case the pipeline crashes at `Define_RT_train`, remove the RT prediction folder from `results/` and re-execute the rule. Most likely, it will work now.

### Check job status
Enter `squeue` on any of the nodes except the Mascot server to check which jobs are currently running.
You can track the status of your job via the outfile (`data/SPIsnake/outfiles/*.out`). Alternatively, you can monitor the node occupancy by `ssh`-ing on the compute nodes and monitoring their activity.
Finally, when the pipeline finished check if all desired output files are present.

