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
    prot_split = expand(join(dir_DB_exhaustive, ".Split_proteome_chunks_{proteome}.done"), 
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
        directory=dir_DB_exhaustive
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


######################### DEV: ON #########################
# checkpoint code to read command data.frame:

class Checkpoint_MakePattern:
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
        checkpoints.check_Split_proteomes.get(**w)

        # expand pattern
        filename = self.get_filename()

        pattern = expand(self.pattern, filename=filename, **w)
        return pattern

rule make_all_files:
    input:
        Checkpoint_MakePattern(join(dir_DB_exhaustive, ".Generate_{filename}.done"))        
    output:
        join(dir_DB_exhaustive, ".Generate_peptides.done")



rule Generate_peptides:
    input: 
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv")
    output:
        peptide_done = touch(join(dir_DB_exhaustive, ".Generate_{filename}.done"))
    benchmark: 
        join(benchmarks, "Generate_peptides_{filename}.json")
    log: 
        join(logs, "Generate_peptides_{filename}.txt")
    conda: 
        "R_env.yaml"
    resources: # 1 per node
        load = 100 
    params:
        directory=dir_DB_exhaustive,
        Filename = "{filename}"
    script:
        "02_3_Generate_peptides.R"  