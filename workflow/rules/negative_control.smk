

rule fit_group_phi:
    input: 
        data=f"{config['input_dir']}/{{data}}.RData",
        modelfile="../src/models/syn_groups_phi.stan"
    output:
        directory(f"{RES_DIR}/{{data}}/negative_control")
    log:
        "logs/{data}/negative_control.log"
    resources:
        h_rt="2:00:00",
        mem_mb=8000,
        disk_mb=8000
    shell:
        f"Rscript scripts/negative_control.R {{input.data}} {{input.modelfile}} {{output}}"