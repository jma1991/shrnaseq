
rule expression_heatmap:
    input:
        rds=["results/{contrast}-glmLRT.rds", 
        "results/filter_hairpins.rds", 
        "results/corrected_counts.rds"]
    output:
        plot=["plots/{contrast}-expression-heatmap.png",
        "plots/{contrast}-corrected-expression-heatmap.png"]
    params:
        FC=config["FC"]
    log:
        out = "logs/{contrast}-expression-heatmap.out",
        err = "logs/{contrast}-expression-heatmap.err" 
    message:
        "Visualise differential expression across groups"
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/expression_heatmap.R"

rule volcano_plot:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        plot="plots/{contrast}-volcano-plot.png"
    log:
        out = "logs/{contrast}-volcano-plot.out",
        err = "logs/{contrast}-volcano-plot.err" 
    message:
        "Generate volcano plot to visualise relationship between magnitude and strength of evidence"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/volcano_plot.R"
   
rule hairpin_histogram:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        plot="plots/{contrast}-hairpin-histogram.png"
    log:
        out = "logs/{contrast}-hairpin-histogram.out",
        err = "logs/{contrast}-hairpin-histogram.err"
    message:
        "Visualise distribution of hairpin p values"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/hairpin_histogram.R"

rule plotSMEAR:
    input:
        rds=["results/{contrast}-glmLRT.rds", "results/{contrast}-FDR_hairpins.rds"]
    output:
        plot="plots/{contrast}-plotSmear.png"
    params:
        FC=config["FC"]
    log:
        out = "logs/{contrast}-plotSmear.out",
        err = "logs/{contrast}-plotSmear.err"
    message:
        "Visualise logFC against logCPM with top FDR hairpins highlighted"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/plotSMEAR.R"       
