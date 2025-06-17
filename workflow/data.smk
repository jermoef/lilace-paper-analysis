# RES_DIR = f"{config['res_dir']}/selfphi_{config['selfphi']}"

RES_DIR = f"{config['res_dir']}"
SCRATCH_DIR = f"{config['scratch_dir']}"
# include: "rules/negative_control.smk"
include: "rules/run.smk"
include: "rules/data_comparison.smk"

desired_files = []
for data, baselines in zip(config['datasets'], config['baselines']):
    for baseline in baselines:
        # data fit
        desired_files.append(f"{RES_DIR}/{data}/baseline_{baseline}/fitted_df.RData")
        # masked fit
        for p_masked in config['prop_masked']:
            for i in range(0, config['n_masked_iter']):
                masked_suffix = f"masked_{p_masked}_{i}"
                desired_files.append(f"{RES_DIR}/{data}_{masked_suffix}/baseline_{baseline}_{masked_suffix}/fitted_df.RData")
        # FDR stats
        # desired_files.append(f"{RES_DIR}/{data}/baseline_{baseline}/FDR.tsv")
        # also desire plots (use last plot to be made as desired output)
        for model in config['models']:
            desired_files.append(f"{RES_DIR}/{data}/baseline_{baseline}/{model}_plots/variant_props_by_score.html")

rule all:
    input: desired_files