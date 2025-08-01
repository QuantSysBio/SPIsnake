"""Snakemake rules for the workflow"""

from glob import glob
from os.path import join
import random
import string

seed = 1
random.seed(seed)

# Wildcards: proteome pre-processing
Master_table = (
    pd.read_csv(features["Master_table"], sep=",")
)

def get_cluster_proteomes_input(Master_table) :
    prot_cluster = expand(join(dir_cluster, "{proteome}/{proteome}_cluster.tsv"), 
            proteome = Master_table["Proteome"].unique())
    return(prot_cluster)


def get_split_proteomes_input(Master_table) :
    prot_split = expand(join(dir_DB_Fasta_chunks, "{proteome}.parquet"), 
            proteome = Master_table["Proteome"].unique())
    return(prot_split)


rule Index_proteome:
    """
    Use linclust to order proteins by sequence similarity
    """
    input:
        proteome = join(dir_reference, '{proteome}.fasta')
    output:
        prot_index = join(dir_reference, '{proteome}.fasta.fai')
    benchmark: 
        join(benchmarks, "Index_proteome_{proteome}.json")
    log: 
        join(logs, "Index_proteome_{proteome}.txt")
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_cluster_proteome=join(dir_cluster, "{proteome}/{proteome}"),
        dir_tmp=join(dir_cluster, "{proteome}/tmp")
    shell:
        "sed 's, ,_,g' -i {input.proteome} ;\
        samtools faidx {input.proteome} &> {log}"


rule Cluster_proteome:
    """
    Use linclust to order proteins by sequence similarity
    """
    input:
        proteome = join(dir_reference, '{proteome}.fasta')
    output:
        prot_cluster = join(dir_cluster, "{proteome}/{proteome}_cluster.tsv")
    benchmark: 
        join(benchmarks, "Cluster_proteome_{proteome}.json")
    log: 
        join(logs, "Cluster_proteome_{proteome}.txt")
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_cluster_proteome=join(dir_cluster, "{proteome}/{proteome}"),
        dir_tmp=join(dir_cluster, "{proteome}/tmp")
    shell:
        "mmseqs easy-linclust -e 1.000E-03 \
        --spaced-kmer-mode 1 --spaced-kmer-pattern 110101 \
        {input.proteome} {params.dir_cluster_proteome} {params.dir_tmp} \
         &> {log}"


rule Convert_RAW:
    """
    Convert RAW files to MGF and Extract MW & RT information 
    Extract RT and peptide sequences for RT predictor training
    """
    input:
        Experiment_design = features["Experiment_design"],
        PEAKS_PSM = features["PEAKS_PSM"],
        allowed_modifications = features["Experiment_design"]
    output:
        Convert_RAW_done = touch(join(dir_DB_exhaustive, ".Convert_RAW.done"))
    benchmark: 
        join(benchmarks, "Convert_RAW.json")
    log: 
        join(logs, "Convert_RAW.txt")
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        ThermoRawFileParser = features["RT_filter"]["ThermoRawFileParser"],
        dir_mass_list = dir_mass_list,
        dir_RAW = dir_RAW,
        dir_RT_calibration = dir_RT_calibration,
        dir_DB_exhaustive = dir_DB_exhaustive,
        max_cpus = config["max_cpus"]
    script:
        "01_01_Convert_RAW.R"
        

rule Split_proteome_chunks:
    """
    Splits input proteome into chunks of approx similar volume for peptide generation
    """
    input:
        proteome = join(dir_reference, '{proteome}.fasta'),
        prot_index = join(dir_reference, '{proteome}.fasta.fai'),
        prot_cluster = join(dir_cluster, "{proteome}/{proteome}_cluster.tsv"),
        Master_table = features["Master_table"]
    output:
        Split_prot_cluster = join(dir_DB_Fasta_chunks, "{proteome}.parquet")
    benchmark: 
        join(benchmarks, "Split_proteome_chunks_{proteome}.json")
    log: 
        join(logs, "Split_proteome_chunks_{proteome}.txt")
    resources:
        ncpus = 1,
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        min_protein_length = features["DB"]["min_protein_length"],
        max_protein_length = features["DB"]["max_protein_length"],
        replace_I_with_L = features["DB"]["replace_I_with_L"],
        directory=dir_DB_Fasta_chunks
    script:
        "02_1_Split_proteome_chunks.R"


rule Expand_Master_table:
    """
    Create a table to control PCP/PSP generation
    """
    input: 
        Split_chunks = get_split_proteomes_input(Master_table),
        Master_table = features["Master_table"],
        Convert_RAW_done = join(dir_DB_exhaustive, ".Convert_RAW.done")
    output:
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv")
    benchmark: 
        join(benchmarks, "Expand_Master_table.json")
    log: 
        join(logs, "Expand_Master_table.txt")
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        directory=dir_DB_exhaustive,
        dir_RT_calibration=dir_RT_calibration,
        max_cpus = config["max_cpus"]
    script:
        "02_2_Expand_Master_table.R"        


# Checkpoint for rules that depend on proteome chunks
# Expand_Master_table.R makes a table with single row per combination
# of proteome chunk and peptide generation parameters as definded by Master_table.csv
checkpoint check_Master_table_expanded:
    input:
        join(dir_DB_exhaustive, "Master_table_expanded.csv")
    output:
        touch(join(dir_DB_exhaustive, ".Expand_Master_table.done"))


# ---------------------------------------- Generate peptides ------------------------------------
# checkpoint code to read command data.frame:
class Checkpoint_Master_table_expanded:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(Master_table_expanded) :
        Master_table_expanded = pd.read_csv(join(dir_DB_exhaustive, "Master_table_expanded.csv"), sep=",")
        filename = Master_table_expanded["filename"]
        return(filename)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Expand_Master_table'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_Master_table_expanded.get(**w)

        # expand pattern
        filename = self.get_filename()

        pattern = expand(self.pattern, filename=filename, **w)
        return pattern


rule make_all_files:
    input:
        Checkpoint_Master_table_expanded(join(dir_DB_exhaustive, ".{filename}.done"))
    output:
        touch(join(dir_DB_exhaustive, ".Generate_peptides.done"))
    resources:
        ncpus = 1,
        mem = config["min_mem"]


rule Generate_peptides:
    input: 
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        RT_Performance_df = join(dir_RT_prediction, "RT_Performance.csv"),
        Experiment_design = features["Experiment_design"]
    output:
        Generate_peptides_done = touch(join(dir_DB_exhaustive, ".{filename}.done"))
    benchmark: 
        join(benchmarks, "Generate_peptides_{filename}.json")
    log: 
        join(logs, "Generate_peptides_{filename}.txt")
    resources:
        ncpus = config["cpus_critical"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_Fasta_chunks=dir_DB_Fasta_chunks,
        dir_PSP_indices=dir_PSP_indices,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_DB_out=dir_DB_out,
        min_protein_length=features["DB"]["min_protein_length"],
        max_protein_length=features["DB"]["max_protein_length"],
        replace_I_with_L=features["DB"]["replace_I_with_L"],
        max_variable_PTM=features["DB"]["max_variable_PTM"],
        generate_spliced_PTMs=features["DB"]["generate_spliced_PTMs"],
        PTM_chunk=features["DB"]["PTM_chunk"],
        netMHCpan_path=features["DB"]["netMHCpan_path"],
        netMHCpan_chunk=features["DB"]["netMHCpan_chunk"],
        AA_index_length=features["DB"]["AA_index_length"],
        compress_CSV=features["FASTA_outputs"]["compress_CSV"],
        FASTA_outputs_unfiltered=features["FASTA_outputs"]["unfiltered"],
        FASTA_outputs_MW_filtered=features["FASTA_outputs"]["MW_filtered"],
        FASTA_outputs_MW_filtered_PTM=features["FASTA_outputs"]["MW_filtered_PTM"],
        FASTA_outputs_MW_RT_filtered=features["FASTA_outputs"]["MW_RT_filtered"],
        FASTA_outputs_MW_RT_IC50_filtered=features["FASTA_outputs"]["MW_RT_IC50_filtered"],
        method=features["RT_filter"]["method"],
        max_cpus = config["max_cpus"]
    script:
        "02_4_SPIsnake_main.R"


rule Define_arrow_aggregation:
    input: 
        Generate_peptides = join(dir_DB_exhaustive, ".Generate_peptides.done")
    output:
        # DB_pep_map = join(dir_DB_out, "DB_pep_map.csv"),
        # DB_PTM_mz = join(dir_DB_out, "DB_PTM_mz.csv"),
        arrow_batch_definition = join(dir_DB_out, "arrow_batch_definition.csv"),
        arrow_batch_definition_map = join(dir_DB_out, "arrow_batch_definition_map.csv")
    benchmark: 
        join(benchmarks, "Define_arrow_aggregation.json")
    log: 
        join(logs, "Define_arrow_aggregation.txt")
    resources:
        ncpus = config["cpus_critical"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_DB_exhaustive = dir_DB_exhaustive,
        dir_DB_PTM_mz = dir_DB_PTM_mz,
        dir_DB_out = dir_DB_out,
        duckdb_max_filesize = features["DB"]["duckdb_max_filesize"],
        max_cpus = config["max_cpus"]
    script:
        "02_4_3_Define_arrow_aggregation.R"


# Checkpoint for arrow aggregation
checkpoint check_arrow_batch_definition:
    input:
        arrow_batch_definition = join(dir_DB_out, "arrow_batch_definition.csv"),
        arrow_batch_definition_map = join(dir_DB_out, "arrow_batch_definition_map.csv")
    output:
        arrow_batch_definition_done = touch(join(dir_DB_out, ".arrow_batch_definition.done"))


# checkpoint code to read command data.frame:
class Checkpoint_arrow_map_aggregation:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(arrow_batch_table) :
        arrow_batch_table = pd.read_csv(join(dir_DB_out, "arrow_batch_definition_map.csv"), sep=',')
        aggregation_batch = arrow_batch_table["aggregation_batch"].drop_duplicates()
        return(aggregation_batch)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Define_arrow_aggregation'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_arrow_batch_definition.get(**w)

        # expand pattern
        aggregation_batch = self.get_filename()

        pattern = expand(self.pattern, aggregation_batch=aggregation_batch, **w)
        return pattern


# checkpoint code to read command data.frame:
class Checkpoint_arrow_aggregation:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(arrow_batch_table) :
        arrow_batch_table = pd.read_csv(join(dir_DB_out, "arrow_batch_definition.csv"), sep=',')
        aggregation_batch = arrow_batch_table["aggregation_batch"].drop_duplicates()
        return(aggregation_batch)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Define_arrow_aggregation'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_arrow_batch_definition.get(**w)

        # expand pattern
        aggregation_batch = self.get_filename()

        pattern = expand(self.pattern, aggregation_batch=aggregation_batch, **w)
        return pattern


rule Aggregare_arrow:
    input:
        arrow_batch_definition = join(dir_DB_out, "arrow_batch_definition.csv")
    output:
        arrow_batch_done = touch(join(dir_DB_out, ".peptides_{aggregation_batch}.done"))
    benchmark: 
        join(benchmarks, "Aggregare_arrow_{aggregation_batch}.json")
    log: 
        join(logs, "Aggregare_arrow_{aggregation_batch}.txt")
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_DB_out=dir_DB_out,
        duckdb_max_filesize = features["DB"]["duckdb_max_filesize"],
        duckdb_max_retries = features["DB"]["duckdb_max_retries"],
        duckdb_RAM = features["DB"]["duckdb_RAM"],
        max_cpus = config["max_cpus"]
    script:
        "02_4_4_Aggregare_arrow.R"


rule Aggregare_arrow_peptides:
    input:
        Checkpoint_arrow_aggregation(join(dir_DB_out, ".peptides_{aggregation_batch}.done"))
    output:
        touch(join(dir_DB_out, ".Aggregare_arrow_peptides.done"))
    resources:
        ncpus = 1,
        mem = config["min_mem"]
        
        
rule Aggregare_peptide_mapping:
    input:
        # DB_pep_map = join(dir_DB_out, "DB_pep_map.csv"),
        arrow_batch_definition_map = join(dir_DB_out, "arrow_batch_definition_map.csv")
    output:
        arrow_map_batch_done = touch(join(dir_DB_out, ".mapping_{aggregation_batch}.done"))
    benchmark: 
        join(benchmarks, "Aggregare_arrow_map_{aggregation_batch}.json")
    log: 
        join(logs, "Aggregare_arrow_map_{aggregation_batch}.txt")
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_DB_out=dir_DB_out,
        duckdb_max_filesize = features["DB"]["duckdb_max_filesize"],
        duckdb_max_retries = features["DB"]["duckdb_max_retries"],
        duckdb_RAM = features["DB"]["duckdb_RAM"],
        max_cpus = config["max_cpus"]
    script:
        "02_4_6_Aggregare_arrow_peptide_mapping.R"


rule Aggregare_arrow_mapping:
    input:
        Checkpoint_arrow_aggregation(join(dir_DB_out, ".mapping_{aggregation_batch}.done"))
    output:
        touch(join(dir_DB_out, ".Aggregare_peptide_mapping.done"))
    resources:
        ncpus = 1,
        mem = config["min_mem"]
        

rule Aggregare_FASTA:
    input: 
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        Experiment_design = features["Experiment_design"],
        Aggregare_arrow_done = join(dir_DB_out, ".Aggregare_arrow_peptides.done")
    output:
        touch(join(dir_DB_out, ".Aggregare_FASTA.done"))
    benchmark: 
        join(benchmarks, "Aggregare_FASTA.json")
    log: 
        join(logs, "Aggregare_FASTA.txt")
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_DB_out=dir_DB_out,
        duckdb_max_filesize = features["DB"]["duckdb_max_filesize"],
        duckdb_max_retries = features["DB"]["duckdb_max_retries"],
        duckdb_RAM = features["DB"]["duckdb_RAM"],
        FASTA_outputs_unfiltered = features["FASTA_outputs"]["unfiltered"],
        FASTA_outputs_MW_filtered = features["FASTA_outputs"]["MW_filtered"],
        FASTA_outputs_MW_filtered_PTM = features["FASTA_outputs"]["MW_filtered_PTM"],
        FASTA_outputs_MW_RT_filtered = features["FASTA_outputs"]["MW_RT_filtered"],
        FASTA_outputs_MW_RT_IC50_filtered = features["FASTA_outputs"]["MW_RT_IC50_filtered"],
        max_cpus = config["max_cpus"]
    script:
        "02_4_4_Aggregare_FASTA.R"  


rule Aggregare_Stats:
    input: 
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        Experiment_design = features["Experiment_design"],
        Aggregare_arrow_done = join(dir_DB_out, ".Aggregare_arrow_peptides.done")
    output:
        touch(join(dir_DB_out, ".Aggregare_Stats.done"))
    benchmark: 
        join(benchmarks, "Aggregare_Stats.json")
    log: 
        join(logs, "Aggregare_Stats.txt")
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"],
        time = config["max_time"]
    params:
        max_cpus = config["max_cpus"],
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_DB_out=dir_DB_out,
        strata_sizes=features["Statistics"]["strata_sizes"],
        filtering_sizes=features["Statistics"]["filtering_sizes"]
    script:
        "02_4_5_Statistics.R"
