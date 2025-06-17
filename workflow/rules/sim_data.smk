

rule get_dataset_params:
    input:
        data=f"{config['sim_input_dir']}/{{data}}.RData",
        desc=f"{config['sim_input_dir']}/{{data}}.yaml"
    output:
        param_file=f"{SIM_RES_DIR}/param_data/param_estimated_{{data}}.RData",
        param_yaml=f"{SIM_RES_DIR}/param_data/param_estimated_{{data}}.yaml"
    resources:
        h_rt="4:00:00",
        mem_mb=8000,
        disk_mb=8000
    shell:
        f"Rscript scripts/get_sim_params.R {{input.data}} {{input.desc}} {{output.param_file}} {{output.param_yaml}}"

rule generate_sim_data:
    input:
        param_file=f"{SIM_RES_DIR}/param_data/param_estimated_{{data}}.RData",
        param_yaml=f"{SIM_RES_DIR}/param_data/param_estimated_{{data}}.yaml"
    output:
        obs_df=f"{SIM_RES_DIR}/sim_data/sim_{{data}}_{{n_pos}}pos_{{param}}_{{param_val}}_{{iter}}.RData",
        obs_df_rescaled=f"{SIM_RES_DIR}/sim_data/sim_{{data}}_{{n_pos}}pos_{{param}}_{{param_val}}_{{iter}}_rescaled.RData",
        sim_obj=f"{SIM_RES_DIR}/param_data/sim_obj_{{data}}_{{n_pos}}pos_{{param}}_{{param_val}}_{{iter}}.RData",
        yaml=f"{SIM_RES_DIR}/sim_data/sim_{{data}}_{{n_pos}}pos_{{param}}_{{param_val}}_{{iter}}.yaml",
        yaml_rescaled=f"{SIM_RES_DIR}/sim_data/sim_{{data}}_{{n_pos}}pos_{{param}}_{{param_val}}_{{iter}}_rescaled.yaml",
    resources:
        h_rt="8:00:00",
        mem_mb=8000,
        disk_mb=8000
    shell:
        f"Rscript ../src/sim_param.R {{input.param_file}} {{input.param_yaml}} {{wildcards.n_pos}} {{wildcards.param}} {{wildcards.param_val}} {{output.obs_df}} {{output.obs_df_rescaled}} {{output.sim_obj}} {{wildcards.iter}}"


# # plot param dist (effect sizes colored by type + text proportion) and simulated data
# rule plot_sim_data:
#     ...


