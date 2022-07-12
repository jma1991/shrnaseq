rule processAmplicons:
    input: 
        index=config["samples"],
        hairpins=config["hairpins"],
        fastq=config["fastq"]
    output: 
        rds="results/processAmplicons.rds"
    message:
        "Process raw read data from genetic sequencing screens"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/processAmplicons.R"

rule filter_hairpins:
    input:
        rds="results/processAmplicons.rds"
    output:
        rds="results/filter_hairpins.rds",
        plot="plots/counts-index-hairpins.png"
    message:
        "Filter hairpins with counts of > 0.5 counts per million in at least 3 samples, and plots of counts per index and hairpins"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/filter_hairpins.R"

rule MDS_plot:
    input:
        rds="results/filter_hairpins.rds"
    output:
        plot="plots/MDS-plot.png"
    message:
        "Multidimensional Scaling plot to visualise relationship between samples"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/MDS_plot.R"

rule model_matrix:
    input:
        rds="results/filter_hairpins.rds"
    output:
        rds="results/model_matrix.rds"
    message:
        "Generate model matrix for edgeR analysis"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/model_matrix.R"

rule contrasts_matrix:
    input:
        rds="results/model_matrix.rds"
    output:
        rds="results/contrasts_matrix.rds"
    message:
        "Generate matrix for group contrasts"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/contrasts_matrix.R"

rule diff_rep_analysis:
    input:
        rds=["results/filter_hairpins.rds", "results/model_matrix.rds"]
    output:
        rds="results/diff_rep_analysis.rds"
    message: 
        "Perform differential representation analysis"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/diff_rep_analysis.R" 

rule BVC_plot:
    input:
        rds="results/diff_rep_analysis.rds"
    output:
        plot="plots/BCV-plot.png"
    message:
        "Visualise Biological Coefficient of Variation against read abundance"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/BVC_plot.R" 

rule PCA_plot:
    input:
        rds="results/filter_hairpins.rds"
    output:
        plot="plots/PCA-plot.png"
    message:
        "Visualise relationships between first 2 principal components"
    conda:
        "../envs/pca.yaml"
    script:
        "../scripts/PCA_plot.R"
    

rule sample_dist_heatmap:
    input:
        rds="results/filter_hairpins.rds"
    output:
        plot="plots/sample-dist-heatmap.png"
    message:
        "Heatmap of sample distances"
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/sample_dist_heatmap.R"
    
rule glmFit:
    input:
        rds=["results/model_matrix.rds","results/diff_rep_analysis.rds"]
    output:
        rds="results/glmFit.rds"
    message:
        "Fit negative bionomial GLM"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmFit.R" 

rule glmLRT:
    input:
        rds=["results/contrasts_matrix.rds", "results/glmFit.rds"]
    output:
        rds="results/glmLRT.rds"
    params:
        contrast=config["contrast"]
    message:
        "Perform likelihood ratio test on GLM"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmLRT.R" 

rule expression_heatmap:
    input:
        rds=["results/filter_hairpins.rds","results/glmLRT.rds"]
    output:
        plot="plots/expression-heatmap.png"
    params:
        FC=config["FC"]
    message:
        "Visualise differential expression across groups"
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/expression_heatmap.R"

rule volcano_plot:
    input:
        rds="results/glmLRT.rds"
    output:
        plot="plots/volcano-plot.png"
    message:
        "Generate volcano plot to visualise relationship between magnitude and strength of evidence"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/volcano_plot.R"
    
rule hairpin_histogram:
    input:
        rds="results/glmLRT.rds"
    output:
        plot="plots/hairpin-histogram-plot.png"
    message:
        "Visualise distribution of hairpin p values"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/hairpin_histogram.R"

rule top_hairpins:
    input:
        rds="results/glmLRT.rds"
    output:
        txt="results/top-ranked-hairpins.txt"
    message:
        "Generate table of the top differentially expressed hairpins"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/top_hairpins.R"

rule FDR_hairpins:
    input:
        rds="results/glmLRT.rds"
    output:
        txt="results/FDR-sig-hairpins.txt",
        rds="results/FDR_hairpins.rds"
    params:
        threshold=config["FDR"]
    message:
        "Highlight and generate table of hairpins with FDR < 0.05"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/FDR_hairpins.R"

rule plotSMEAR:
    input:
        rds=["results/glmLRT.rds", "results/FDR_hairpins.rds"]
    output:
        plot="plots/plotSmear.png"
    message:
        "Visualise logFC against logCPM with top FDR hairpins highlighted"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/plotSMEAR.R"       


