
rule make_masked:
    input:
        data=f"{config['input_dir']}/{{data}}.RData",
        desc=f"{config['input_dir']}/{{data}}.yaml"
    output:
        masked=f"{config['input_dir']}/{{data}}_masked_{{prop_masked}}_{{iter}}.RData",
        desc=f"{config['input_dir']}/{{data}}_masked_{{prop_masked}}_{{iter}}.yaml"
    shell:
        f"Rscript scripts/make_masked.R {{input.data}} {{input.desc}} {{output.masked}} {{output.desc}} {{wildcards.prop_masked}} {{wildcards.iter}}"


# run each stan model specified as separate jobs
rule run_model:
    input:
        data=f"{config['input_dir']}/{{data}}.RData",
        baseline=f"{config['input_dir']}/{{baseline}}.RData",
        modelfile="../src/models/{model}.stan",
        baseline_modelfile="../src/models/baseline.stan",
        desc=f"{config['input_dir']}/{{data}}.yaml"
    output: 
        stanfit=f"{SCRATCH_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_stanfit.RData",
        out_df=f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_df.RData",
        input_file=f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_input_df.RData"
    params:
        plot_dir=directory(f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_plots")
    log:
        "logs/{data}/baseline_{baseline}/{model}.log"
    resources:
        h_rt="23:30:00",
        mem_mb=8000,
        disk_mb=8000,
        # highp=""
    # resources:
    #     h_rt="45:30:00",
    #     mem_mb=16000,
    #     disk_mb=16000,
    #     highp=""
    shell:
        f"Rscript scripts/run_model.R {{input.modelfile}} {{input.baseline_modelfile}} {{input.data}} {{input.baseline}} {{input.desc}} {config['selfphi']} {{output.stanfit}} {{output.out_df}} {{output.input_file}} {{params.plot_dir}}"

rule summarize_model:
    input:
        modelfile="../src/models/{model}.stan",
        stanfit=f"{SCRATCH_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_stanfit.RData",
        out_df=f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_df.RData"
    output:
        f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_plots/variant_props_by_score.html" # last plot/file to be generated
    params:
        plot_dir=directory(f"{RES_DIR}/{{data}}/baseline_{{baseline}}/{{model}}_plots")
    # resources:
    #     h_rt="6:00:00",
    #     mem_mb=48000,
    #     disk_mb=32000
    resources:
        h_rt="12:30:00",
        mem_mb=8000,
        disk_mb=8000,
    shell:
        f"Rscript scripts/summarize_model.R {{input.modelfile}} {{input.stanfit}} {{input.out_df}} {{params.plot_dir}}"

# run all other methods specified in config in one job bc they prob won't take as long
rule run_methods:
    input:
        data=f"{config['input_dir']}/{{data}}.RData",
        desc=f"{config['input_dir']}/{{data}}.yaml",
    output:
        f"{RES_DIR}/{{data}}/methods_df.RData"
    log:
        "logs/{data}/methods.log"
    resources:
        h_rt="12:00:00",
        mem_mb=8000,
        disk_mb=8000
    shell:
        f"Rscript scripts/run_methods.R {{input.data}} {{input.desc}} {','.join(config['methods'])} {RES_DIR}/{{wildcards.data}} {{output}} \
        {config['SHELL']} {config['RCFILE']}"


rule combine_fits:
    input:
       expand("%s/{{data}}/baseline_{{baseline}}/{model}_df.RData" % RES_DIR, model=config['models']),
       f"{RES_DIR}/{{data}}/methods_df.RData"
    output:
        f"{RES_DIR}/{{data}}/baseline_{{baseline}}/fitted_df.RData",
    log:
        "logs/{data}/baseline_{baseline}/combine.log"
    resources:
        h_rt="0:10:00",
        mem_mb=8000,
        disk_mb=8000
    shell:
        f"Rscript scripts/combine_fits.R {{input}} {{output}}"
