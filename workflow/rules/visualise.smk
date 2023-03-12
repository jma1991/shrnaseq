# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule plotHeatmap: 
    input:
        rds = "results/calcNormFactors.rds",
        tsv = "results/{contrast}.topTags.tsv"
    output:
        png = "plots/{contrast}.plotHeatmap.png"
    params:
        ntop = 50
    log:
        out = "logs/{contrast}.plotHeatmap.out",
        err = "logs/{contrast}.plotHeatmap.err" 
    message:
        "Plot expression of top-ranked genes for contrast: {wildcards.contrast}"
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/plotHeatmap.R"

rule plotVolcano:
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
   
rule plotPValue:
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

rule plotRank:
    input:
        tsv = "results/{contrast}.combineTests.tsv"
    output:
        png = "results/{contrast}.plotRank.png"
    log:
        out = "logs/{contrast}.plotRank.out",
        err = "logs/{contrast}.plotRank.err"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plotRank.R"
