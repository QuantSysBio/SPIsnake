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
            proteome = Master_table["Proteome"])
    return(prot_cluster)


def get_split_proteomes_input(Master_table) :
    prot_split = expand(join(dir_DB_exhaustive, ".Split_proteome_chunks_{proteome}.done"), 
            proteome = Master_table["Proteome"])
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
        "mmseqs.yaml"
    params:
        dir_cluster_proteome=join(dir_cluster, "{proteome}/{proteome}")
    shell:
        "mkdir tmp & \
        mmseqs easy-linclust -e Inf -c 0.01 \
        --spaced-kmer-mode 1 --spaced-kmer-pattern 110101 \
        {input.proteome} {params.dir_cluster_proteome} tmp \
         &> {log}"



rule Split_proteome_chunks:
    """
    Splits input proteome into chunks of approx similar volume for DB construction
    """
    input:
        proteome = join(dir_reference, '{proteome}.fasta'),
        prot_cluster = join(dir_cluster, "{proteome}/{proteome}_cluster.tsv"),
        Master_table = features["Master_table"],
        functions = "src/snakefiles/functions.R"
    output:
        Split_dir = directory(join(dir_DB_exhaustive, "{proteome}")),
        Split_proteomes = touch(join(dir_DB_exhaustive, ".Split_proteome_chunks_{proteome}.done"))
    benchmark: 
        join(benchmarks, "Split_proteome_chunks_{proteome}.json")
    log: 
        join(logs, "Split_proteome_chunks_{proteome}.txt")
    conda: 
        "R_env.yaml"
    params:
        n=config["max_cores"],
        min_protein_length = features["DB"]["min_protein_length"],
        max_protein_length = features["DB"]["max_protein_length"],
        maxE = features["DB"]["maxE"],
        dir_DB_exhaustive=dir_DB_exhaustive
    script:
        "02_1_Split_proteome_chunks.R"


rule Expand_Master_table:
    input: 
        Split_chunks = get_split_proteomes_input(Master_table),
        Master_table = features["Master_table"]
    output:
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv")
    benchmark: 
        join(benchmarks, "Expand_Master_table.json")
    log: 
        join(logs, "Expand_Master_table.txt")
    conda: 
        "R_env.yaml"
    params:
        dir_DB_exhaustive=dir_DB_exhaustive
    script:
        "02_2_Expand_Master_table.R"        


# Checkpoint for rules that depend on proteome chunks"
# Expand_Master_table.R makes a table with single row per combination
# of proteome chunk and peptide generation parameters as definded by Master_table.csv
checkpoint check_Split_proteomes:
    input:
        join(dir_DB_exhaustive, "Master_table_expanded.csv")
    output:
        touch(join(dir_DB_exhaustive, ".Expand_Master_table.done"))


