rule clone_AutoRT:
    output:
        AutoRT = "bin/AutoRT/autort.py"
    benchmark: 
        join(benchmarks, "clone_AutoRT.json")
    log: 
        join(logs, "clone_AutoRT.txt")
    conda: 
        "git.yaml"
    params:
        url = "https://github.com/bzhanglab/AutoRT"
    shell:
        "git clone {params.url} bin/AutoRT"


rule Define_RT_train:
    input:
        AutoRT = "bin/AutoRT/autort.py",
        Experiment_design = features["Experiment_design"]
    output:
        cmd_RT_train = join(dir_RT_prediction, "cmd_RT_train.csv"),
        cmd_RT_test = join(dir_RT_prediction, "cmd_RT_test.csv")
    benchmark: 
        join(benchmarks, "Define_RT_train.json")
    log: 
        join(logs, "Define_RT_train.txt")
    conda: 
        "R_env_reticulate.yaml"
    params:
        method=features["RT_filter"]["method"],
        n_folds=features["RT_filter"]["n_folds"],
        proportion_train=features["RT_filter"]["proportion_train"],
        dir_RT_calibration=dir_RT_calibration,
        dir_RT_prediction=dir_RT_prediction
    script:
        "04_1_Define_RT_train.R"


# Checkpoint for aggregating generated peptides across proteome chunks
checkpoint check_RT_train:
    input:
        cmd_RT_train = join(dir_RT_prediction, "cmd_RT_train.csv")
    output:
        cmd_RT_train_done = touch(join(dir_RT_prediction, ".cmd_RT_train.done"))


# checkpoint code to read command data.frame:
class Checkpoint_RT_train:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_filename(RT_calibration_table) :
        RT_calibration_table = pd.read_csv(join(dir_RT_prediction, "cmd_RT_train.csv"), sep=',')
        RT_dataset = RT_calibration_table["RT_dataset"]
        return(RT_dataset)

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'Define_RT_train'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_RT_train.get(**w)

        # expand pattern
        RT_dataset = self.get_filename()

        pattern = expand(self.pattern, RT_dataset=RT_dataset, **w)
        return pattern


rule aggregate_RT:
    input:
        checkpoint = Checkpoint_RT_train(join(dir_RT_prediction, "RT_models/{RT_dataset}/model.json")),
        cmd_RT_test = join(dir_RT_prediction, "cmd_RT_test.csv")
    output:
        RT_Performance_df = join(dir_RT_prediction, "RT_Performance.csv")
    benchmark: 
        join(benchmarks, "aggregate_RT.json")
    log: 
        join(logs, "aggregate_RT.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        load = 100 
    params:
        method=features["RT_filter"]["method"]
    script:
        "04_3_aggregate_RT.R"


rule RT_train:
    input: 
        cmd_RT_train = join(dir_RT_prediction, "cmd_RT_train.csv")
    output:
        RT_best_model = join(dir_RT_prediction, "RT_models/{RT_dataset}/model.json")
    benchmark: 
        join(benchmarks, "RT_train_{RT_dataset}.json")
    log: 
        join(logs, "RT_train_{RT_dataset}.txt")
    conda: 
        "R_env_reticulate.yaml"
    resources: # 1 per node at the time
        load = 100 
    script:
        "04_2_RT_train.R"


