rule Define_peptide_aggregation:
    input: 
        Generate_peptides_done = join(dir_DB_exhaustive, ".Generate_peptides.done"),
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv")
    output:
        Peptide_aggregation_table = join(dir_DB_PTM_mz, "Peptide_aggregation_table.csv")
    benchmark: 
        join(benchmarks, "Define_peptide_processing.json")
    log: 
        join(logs, "Define_peptide_processing.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        AA_index_length=features["DB"]["AA_index_length"],
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz
    script:
        "03_1_Define_peptide_aggregation.R"


# Checkpoint for aggregating generated peptides across proteome chunks
checkpoint check_Peptide_aggregation:
    input:
        join(dir_DB_PTM_mz, "Peptide_aggregation_table.csv")
    output:
        Peptide_aggregation_table_done = touch(join(dir_DB_PTM_mz, ".Peptide_aggregation_table.done"))


# checkpoint code to read command data.frame:
class Checkpoint_Peptide_aggregation:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(Peptide_aggregation_table) :
        Peptide_aggregation_table = pd.read_csv(join(dir_DB_PTM_mz, "Peptide_aggregation_table.csv"), sep=",")
        AA_length = Peptide_aggregation_table["AA_length"]
        return(AA_length)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Expand_Master_table'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_Peptide_aggregation.get(**w)

        # expand pattern
        AA_length = self.get_filename()

        pattern = expand(self.pattern, AA_length=AA_length, **w)
        return pattern


rule aggregate_chunks:
    input:
        Checkpoint_Peptide_aggregation(join(dir_DB_PTM_mz, "chunk_aggregation_status/{AA_length}.csv"))        
    output:
        touch(join(dir_DB_PTM_mz, ".Aggregate_peptides.done"))
    resources:
        ncpus = 1,
        mem = config["min_mem"]


rule PTM_mz_RT_matching:
    input: 
        Peptide_aggregation_table = join(dir_DB_PTM_mz, "Peptide_aggregation_table.csv"),
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        Experiment_design = features["Experiment_design"],
        RT_Performance_df = join(dir_RT_prediction, "RT_Performance.csv")
    output:
        chunk_aggregation_status = join(dir_DB_PTM_mz, "chunk_aggregation_status/{AA_length}.csv")
    benchmark: 
        join(benchmarks, "PTM_mz_matching_{AA_length}.json")
    log: 
        join(logs, "PTM_mz_matching_{AA_length}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = config["cpus_critical"],
        mem = config["max_mem"] 
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        method=features["RT_filter"]["method"],
        netMHCpan_chunk=features["DB"]["netMHCpan_chunk"],
        max_variable_PTM=features["DB"]["max_variable_PTM"],
        generate_spliced_PTMs=features["DB"]["generate_spliced_PTMs"],
        fst_compression = features["DB"]["fst_compression"],
        cpus_for_R = config["max_cpus"]
    script:
        "03_2_PTM_mz_RT_matching.R"
