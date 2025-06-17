
SIM_RES_DIR = f"{config['sim_res_dir']}"
RES_DIR = f"{config['res_dir']}"
SCRATCH_DIR = f"{config['scratch_dir']}"


include: "rules/sim_data.smk"
include: "rules/run.smk"

desired_files = []
for data in config['datasets']:
    for pos in config['n_pos_vals']:
        for param, param_vals in zip(config['params'], config['param_vals']):
            for val in param_vals:
                for iter in range(0, config['n_iter']):
                    dataname = f"{data}_{pos}pos_{param}_{val}_{iter}"
                    desired_files.append(f"{SIM_RES_DIR}/sim_data/sim_{dataname}_rescaled.RData")
                    desired_files.append(f"{SIM_RES_DIR}/param_data/sim_obj_{dataname}.RData")
                    desired_files.append(f"{RES_DIR}/sim_{dataname}_rescaled/baseline_sim_{dataname}_rescaled/fitted_df.RData")

rule all:
    input: desired_files