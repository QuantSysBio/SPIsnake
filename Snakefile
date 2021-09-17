shell.prefix("set -euo pipefail;")
shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

container: "docker://snakemake/snakemake:v6.8.0"
min_version("6.0")

config = yaml.load(open("config.yml", "r"), Loader=yaml.FullLoader)
features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/snakefiles/"
include: snakefiles + "folders.py"
include: snakefiles + "reads_wildcards.py"
include: snakefiles + "pre_process.py"
include: snakefiles + "STAR.py"
include: snakefiles + "filter_BAM.py"
include: snakefiles + "stringtie.py"
include: snakefiles + "salmon.py"
include: snakefiles + "transcriptome_2_proteome.py"
include: snakefiles + "multiqc.py"


rule all:
    input:
        ### Trimmed reads
        expand("results/trimmed/{sample}_1.fq.gz", sample = SAMPLES),
        expand("results/trimmed/{sample}_2.fq.gz", sample = SAMPLES),

        ### De-novo transcriptome assembly and protein prediction
        gff_compare="results/stringtie_assemble/gffcompare/gffcmp.combined.gtf",
        proteome_ref_denovo_expressed="results/tr_2_prot/proteome_ref_denovo_expressed.fasta",        

        ### Deep learning coding status prediction
#        classification = "results/RNAsamba_coding_potential/classification.tsv",
#        proteome_ref_denovo_RNAsamba_expressed_nodup="results/tr_2_prot/proteome_ref_denovo_RNAsamba_expressed_nodup.fasta",

        ### Include nORFs
        #expr_prot_nORFs="results/tr_2_prot/proteome_ref_denovo_nORFs_expressed.fasta",
        #expr_prot_nORFs_incl_deep_nodup="results/tr_2_prot/proteome_ref_denovo_RNAsamba_nORFs_expressed_nodup.fasta",
        
        ### Salmon: alignment-based
        #tr_quant_BAM=expand("results/trascriptome_quant/salmon_BAM/{sample}/quant.sf", sample = SAMPLES),
                
        ### Quantification with Stringtie-Ballgown
        # tr_quant_st=expand("results/trascriptome_quant/stringtie/{sample}/{sample}_t_data.ctab", sample=SAMPLES)
        
        ### Salmon: selective alignment
        tr_quant_sa=expand("results/trascriptome_quant/salmon/{sample}/quant.sf", sample = SAMPLES),

        ### MultiQC
        multiqc_report="results/multiqc/multiqc_report.html"
        


### snakemake --dag > dag.dot && dot -Tsvg < dag.dot > dag.svg
### snakemake --filegraph > filegraph.dot && dot -Tsvg < filegraph.dot > filegraph.svg
### snakemake --rulegraph > rulegraph.dot && dot -Tsvg < rulegraph.dot > rulegraph.svg

### snakemake --use-conda --use-singularity -j 1 -r --verbose
### time snakemake --use-singularity --use-conda -j 27 --conda-frontend conda --resources load=100


