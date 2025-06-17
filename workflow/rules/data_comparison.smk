


rule estimate_data_FDR:
    input: 
        f"{RES_DIR}/{{data}}/baseline_{{baseline}}/fitted_df.RData"
    output:
        FDR_out=f"{RES_DIR}/{{data}}/baseline_{{baseline}}/FDR.tsv",
        FDR_calc=f"{RES_DIR}/{{data}}/baseline_{{baseline}}/FDR_calc.tsv",
        FDR_plot=directory(f"{RES_DIR}/{{data}}/baseline_{{baseline}}/FDR_plots")
    log:
        "logs/{data}/baseline_{baseline}/FDR.log"
    resources:
        h_rt="0:30:00",
        mem_mb=4000,
        disk_mb=4000
    shell:
        f"Rscript scripts/FDR_data.R {{input}} {{output.FDR_out}} {{output.FDR_calc}} {{output.FDR_plot}}"