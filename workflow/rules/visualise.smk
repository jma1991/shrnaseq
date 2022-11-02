
rule expression_heatmap: 
    input:
        rds=["results/{contrast}-glmLRT.rds", 
        "results/filter_guideRNAs.rds", 
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
        "../envs/plots.yaml"
    script:
        "../scripts/expression_heatmap.R"

rule volcano_plot:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        plot="plots/{contrast}-volcano-plot.png"
    params:
        FDR=config["FDR"]
    log:
        out = "logs/{contrast}-volcano-plot.out",
        err = "logs/{contrast}-volcano-plot.err" 
    message:
        "Generate volcano plot to visualise relationship between magnitude and strength of evidence"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/volcano_plot.R"
   
rule guideRNA_histogram:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        plot="plots/{contrast}-guideRNA-histogram.png"
    log:
        out = "logs/{contrast}-guideRNA-histogram.out",
        err = "logs/{contrast}-guideRNA-histogram.err"
    message:
        "Visualise distribution of guideRNA p values"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/guideRNA_histogram.R"

rule plotSMEAR:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        plot="plots/{contrast}-plotSmear.png"
    params:
        FC=config["FC"],
        FDR=config["FDR"]
    log:
        out = "logs/{contrast}-plotSmear.out",
        err = "logs/{contrast}-plotSmear.err"
    message:
        "Visualise logFC against logCPM with top FDR guideRNAs highlighted"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/plotSMEAR.R"       
