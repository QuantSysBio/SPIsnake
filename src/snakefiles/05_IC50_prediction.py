rule Define_IC50:
    input:
        Experiment_design = features["Experiment_design"],
        # netMHCpan = "bin/netMHCpan-4.1/netMHCpan",
        mz_RT_aggregation = join(dir_DB_PTM_mz, ".Aggregate_peptides.done")
    output:
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv")
    benchmark: 
        join(benchmarks, "Define_IC50.json")
    log: 
        join(logs, "Define_IC50.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources:
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        dir_IC50=dir_IC50,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        n_netMHCpan_blocks = features["DB"]["n_netMHCpan_blocks"]
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
        cmd_block = netMHCpan_table["cmd_block"].unique()
        return(cmd_block)

    def __call__(self, w):
        global checkpoints

        checkpoints.check_IC50.get(**w)
        cmd_block = self.get_filename()

        pattern = expand(self.pattern, cmd_block=cmd_block, **w)
        return pattern


rule aggregate_IC50_prediction:
    input:
        checkpoint = Checkpoint_IC50(join(dir_IC50, "Seq_stats/{cmd_block}.csv")),
        Experiment_design = features["Experiment_design"],
        Master_table_expanded = join(dir_DB_exhaustive, "Master_table_expanded.csv"),
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv")
    output:
        Summary_stats = join(dir_DB_out, "Stats.csv")
    benchmark: 
        join(benchmarks, "aggregate_IC50.json")
    log: 
        join(logs, "aggregate_IC50.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        dir_DB_exhaustive=dir_DB_exhaustive,
        dir_DB_PTM_mz=dir_DB_PTM_mz,
        dir_IC50=dir_IC50,
        dir_DB_out=dir_DB_out
    script:
        "05_3_aggregate_IC50.R"


rule predict_MHC_affinity:
    input:
        Experiment_design = features["Experiment_design"],
        cmd_netMHCpan = join(dir_IC50, "cmd_netMHCpan.csv")
    output:
        IC50_filter_stats = join(dir_IC50, "Seq_stats/{cmd_block}.csv")
    benchmark: 
        join(benchmarks, "Predict_MHC_affinity_{cmd_block}.json")
    log: 
        join(logs, "Predict_MHC_affinity_{cmd_block}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        ncpus = config["max_cpus"],
        mem = config["max_mem"] 
    params:
        dir_IC50=dir_IC50
    script:
        "05_2_netMHCpan.R"
