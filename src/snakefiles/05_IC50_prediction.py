rule Define_IC50:
    input:
        Experiment_design = features["Experiment_design"],
        netMHCpan = "bin/netMHCpan-4.1/netMHCpan",
        mz_RT_aggregation = join(dir_DB_PTM_mz, ".Aggregate_peptides.done")
    output:
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv"),
        IC50_aggregation_table = join(dir_IC50, "IC50_aggregation_table.csv")
    benchmark: 
        join(benchmarks, "Define_IC50.json")
    log: 
        join(logs, "Define_IC50.txt")
    conda: 
        "R_env.yaml"
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        dir_IC50=dir_IC50,
        dir_DB_PTM_mz=dir_DB_PTM_mz
    script:
        "05_1_Define_IC50.R"


# Checkpoint for predicting IC50 across combinations of peptide chunks and haplotypes
checkpoint check_IC50:
    input:
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv")
    output:
        cmd_netMHCpan_done = touch(join(dir_IC50, ".cmd_netMHCpan.done"))


# checkpoint code to read command data.frame:
class Checkpoint_IC50:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(netMHCpan_table) :
        netMHCpan_table = pd.read_csv(join(dir_IC50, "cmd_netMHCpan.csv"), sep=',')
        Peptide_file = netMHCpan_table["Peptide_file"]
        return(Peptide_file)

    def __call__(self, w):
        global checkpoints

        checkpoints.check_IC50.get(**w)
        Peptide_file = self.get_filename()

        pattern = expand(self.pattern, Peptide_file=Peptide_file, **w)
        return pattern


rule aggregate_IC50_prediction:
    input:
        checkpoint = Checkpoint_IC50(join(dir_IC50, "IC50_filtered_peptides/{Peptide_file}.csv.gz")),
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv")
    output:
        Summary_stats = join(dir_DB_out, "Stats.csv")
    benchmark: 
        join(benchmarks, "aggregate_IC50.json")
    log: 
        join(logs, "aggregate_IC50.txt")
    conda: 
        "R_env.yaml"
    resources: # 1 per node at the time
        load = 100,
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_IC50=dir_IC50,
        dir_DB_out=dir_DB_out
    script:
        "05_3_aggregate_IC50.R"


rule predict_MHC_affinity:
    input: 
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv"),
        netMHCpan = "bin/netMHCpan-4.1/netMHCpan"
    output:
        binders = join(dir_IC50, "IC50_filtered_peptides/{Peptide_file}.csv.gz"),
        IC50_filter_stats = join(dir_IC50, "Seq_stats/{Peptide_file}.csv")
    benchmark: 
        join(benchmarks, "Predict_MHC_affinity_{Peptide_file}.json")
    log: 
        join(logs, "Predict_MHC_affinity_{Peptide_file}.txt")
    conda: 
        "R_env.yaml"
    resources: # 1 per node at the time
        load = 100,
        ncpus = 1,
        mem = 3G 
    container: None
    params:
        dir_IC50=dir_IC50
    script:
        "05_2_netMHCpan.R"
