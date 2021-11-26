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
    prot_split = expand(join(dir_DB_Fasta_chunks, ".Split_proteome_chunks_{proteome}.done"), 
            proteome = Master_table["Proteome"].unique())
    return(prot_split)



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
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"]
    params:
        dir_cluster_proteome=join(dir_cluster, "{proteome}/{proteome}"),
        dir_tmp=join(dir_cluster, "{proteome}/tmp")
    shell:
        "mmseqs easy-linclust -e 1.000E-03 \
        --spaced-kmer-mode 1 --spaced-kmer-pattern 110101 \
        {input.proteome} {params.dir_cluster_proteome} {params.dir_tmp} \
         &> {log}"



rule Split_proteome_chunks:
    """
    Splits input proteome into chunks of approx similar volume for peptide generation
    """
    input:
        proteome = join(dir_reference, '{proteome}.fasta'),
        prot_cluster = join(dir_cluster, "{proteome}/{proteome}_cluster.tsv"),
        Master_table = features["Master_table"],
        functions = "src/snakefiles/functions.R"
    output:
        Split_dir = directory(join(dir_DB_Fasta_chunks, "{proteome}")),
        Split_proteomes = touch(join(dir_DB_Fasta_chunks, ".Split_proteome_chunks_{proteome}.done"))
    benchmark: 
        join(benchmarks, "Split_proteome_chunks_{proteome}.json")
    log: 
        join(logs, "Split_proteome_chunks_{proteome}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = 1,
        mem = config["max_mem"]
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
        Master_table = features["Master_table"]
    output:
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        PSP_indices = join(dir_DB_exhaustive, "PSP_indices.csv")
    benchmark: 
        join(benchmarks, "Expand_Master_table.json")
    log: 
        join(logs, "Expand_Master_table.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = 1,
        mem = config["max_mem"]
    params:
        directory=dir_DB_exhaustive
    script:
        "02_2_Expand_Master_table.R"        


# Checkpoint for rules that depend on proteome chunks
# Expand_Master_table.R makes a table with single row per combination
# of proteome chunk and peptide generation parameters as definded by Master_table.csv
checkpoint check_Split_proteomes:
    input:
        join(dir_DB_exhaustive, "Master_table_expanded.csv")
    output:
        touch(join(dir_DB_exhaustive, ".Expand_Master_table.done"))



# ---------------------------------------- PSP index  ------------------------------------

# checkpoint code to read command data.frame:
class Checkpoint_PSP_indices:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(PSP_indices) :
        PSP_indices = pd.read_csv(join(dir_DB_exhaustive, "PSP_indices.csv"), sep=",")
        PSP_index = PSP_indices["PSP_index"]
        return(PSP_index)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Expand_Master_table'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_Split_proteomes.get(**w)

        # expand pattern
        PSP_index = self.get_filename()

        pattern = expand(self.pattern, PSP_index=PSP_index, **w)
        return pattern


rule Generate_PSP_indices:
    input: 
        PSP_indices = join(dir_DB_exhaustive, "PSP_indices.csv")
    output:
        PSP_index = join(dir_DB_exhaustive, "PSP_indices/{PSP_index}.rds")
    benchmark: 
        join(benchmarks, "Generate_index_{PSP_index}.json")
    log: 
        join(logs, "Generate_index_{PSP_index}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        load = 100,
        ncpus = config["max_cpus"],
        mem = config["max_mem"]
    params:
        directory=dir_DB_exhaustive,
        AA_index_length=features["DB"]["AA_index_length"],
        max_protein_length=features["DB"]["max_protein_length"]
    script:
        "02_3_Generate_PSP_indices.R"  


# ---------------------------------------- Generate peptides ------------------------------------
checkpoint check_Generated_PSP_indices:
    input:
        Checkpoint_PSP_indices(join(dir_DB_exhaustive, "PSP_indices/{PSP_index}.rds"))
    output:
        touch(join(dir_DB_exhaustive, ".Generate_indices.done"))


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
        checkpoints.check_Generated_PSP_indices.get(**w)

        # expand pattern
        filename = self.get_filename()

        pattern = expand(self.pattern, filename=filename, **w)
        return pattern


rule make_all_files:
    input:
        Checkpoint_Master_table_expanded(join(dir_DB_exhaustive, "Seq_stats/{filename}.csv.gz"))
    output:
        touch(join(dir_DB_exhaustive, ".Generate_peptides.done"))
    resources:
        ncpus = 1,
        mem = config["min_mem"]


rule Generate_peptides:
    input: 
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        Generate_indices = join(dir_DB_exhaustive, ".Generate_indices.done")
    output:
        Seq_stats = join(dir_DB_exhaustive, "Seq_stats/{filename}.csv.gz")
    benchmark: 
        join(benchmarks, "Generate_peptides_{filename}.json")
    log: 
        join(logs, "Generate_peptides_{filename}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        load = 100,
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        directory=dir_DB_exhaustive,
        dir_DB_Fasta_chunks=dir_DB_Fasta_chunks,
        AA_index_length=features["DB"]["AA_index_length"],
        max_protein_length=features["DB"]["max_protein_length"]
    script:
        "02_4_Generate_peptides.R"
