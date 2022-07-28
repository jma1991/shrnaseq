rule processAmplicons:
    input: 
        index=config["samples"],
        hairpins=config["hairpins"],
        fastq=config["fastq"]
    output: 
        rds="results/processAmplicons.rds"
    log:
        out = "logs/processAmplicons.out",
        err = "logs/processAmplicons.err"
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
    log:
        out = "logs/filter_hairpins.out",
        err = "logs/filter_hairpins.err"
    message:
        "Filter hairpins with counts of > 0.5 counts per million in at least 3 samples, and plots of counts per index and hairpins"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/filter_hairpins.R"

rule corrected_counts:
    input:
        rds="results/estimateDisp.rds"
    output:
        rds="results/corrected_counts.rds"
    log:
        out = "logs/corrected_counts.out",
        err = "logs/corrected_counts.err"
    message: 
        "Remove batch effect"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/corrected_counts.R" 

rule model_matrix:
    input:
        rds="results/filter_hairpins.rds"
    output:
        rds="results/model_matrix.rds"
    log:
        out = "logs/model_matrix.out",
        err = "logs/model_matrix.err"
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
    log:
        out = "logs/contrasts_matrix.out",
        err = "logs/contrasts_matrix.err"
    message:
        "Generate matrix for group contrasts"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/contrasts_matrix.R"

rule estimateDisp:
    input:
        rds=["results/filter_hairpins.rds", "results/model_matrix.rds"]
    output:
        rds="results/estimateDisp.rds"
    log:
        out = "logs/estimateDisp.out",
        err = "logs/estimateDisp.err"
    message: 
        "Perform differential representation analysis"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/estimateDisp.R" 

rule MDS_plot:
    input:
        rds=["results/filter_hairpins.rds", "results/corrected_counts.rds"]
    output:
        plot=["plots/MDS-plot.png", "plots/corrected-MDS-plot.png"]
    log:
        out = "logs/MDS_plot.out",
        err = "logs/MDS_plot.err"
    message:
        "Multidimensional Scaling plot to visualise relationship between samples"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/MDS_plot.R"

rule PCA_plot:
    input:
        rds=["results/filter_hairpins.rds", "results/corrected_counts.rds"]
    output:
        plot=["plots/PCA-plot.png", "plots/corrected-PCA-plot.png"]
    log:
        out = "logs/PCA_plot.out",
        err = "logs/PCA_plot.err"   
    message:
        "Visualise relationships between first 2 principal components"
    conda:
        "../envs/pca.yaml"
    script:
        "../scripts/PCA_plot.R"

rule sample_dist_heatmap:
    input:
        rds=["results/filter_hairpins.rds",  "results/corrected_counts.rds"]
    output:
        plot=["plots/sample-dist-heatmap.png", "plots/corrected-sample-dist-heatmap.png"]
    log:
        out = "logs/sample_dist_heatmap.out",
        err = "logs/sample_dist_heatmap.err"   
    message:
        "Heatmap of sample distances"
    conda:
        "../envs/heatmap.yaml"
    script:
        "../scripts/sample_dist_heatmap.R"
    
rule BVC_plot:
    input:
        rds="results/estimateDisp.rds"
    output:
        plot="plots/BCV-plot.png"
    log:
        out = "logs/BVC_plot.out",
        err = "logs/BVC_plot.err"    
    message:
        "Visualise Biological Coefficient of Variation against read abundance"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/BVC_plot.R" 

rule glmFit:
    input:
        rds=["results/model_matrix.rds", "results/estimateDisp.rds"]
    output:
        rds="results/glmFit.rds"
    log:
        out = "logs/glmFit.out",
        err = "logs/glmFit.err" 
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
        rds="results/{contrast}-glmLRT.rds"
    params:
        contrast=get_contrast
    log:
        out = "logs/{contrast}-glmLRT.out",
        err = "logs/{contrast}-glmLRT.err" 
    message:
        "Perform likelihood ratio test on GLM"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmLRT.R" 

rule expression_heatmap:
    input:
        rds=["results/{contrast}-glmLRT.rds", 
        "results/filter_hairpins.rds", 
        "results/corrected_counts.rds"]
    output:
        plot=["plots/{contrast}-expression-heatmap.png",
        "plots/corrected-{contrast}-expression-heatmap.png"]
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

rule top_hairpins:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        tsv="results/{contrast}-top-ranked-hairpins.tsv"
    log:
        out = "logs/{contrast}-top-ranked-hairpins.out",
        err = "logs/{contrast}-top-ranked-hairpins.err"
    message:
        "Generate table of the top differentially expressed hairpins"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/top_hairpins.R"

rule FDR_hairpins:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        tsv="results/{contrast}-FDR-sig-hairpins.tsv",
        rds="results/{contrast}-FDR_hairpins.rds"
    params:
        threshold=config["FDR"]
    log:
        out = "logs/{contrast}-FDR-sig-hairpins.out",
        err = "logs/{contrast}-FDR-sig-hairpins.err"
    message:
        "Highlight and generate table of hairpins with FDR < 0.05"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/FDR_hairpins.R"

rule plotSMEAR:
    input:
        rds=["results/{contrast}-glmLRT.rds", "results/{contrast}-FDR_hairpins.rds"]
    output:
        plot="plots/{contrast}-plotSmear.png"
    log:
        out = "logs/{contrast}-plotSmear.out",
        err = "logs/{contrast}-plotSmear.err"
    message:
        "Visualise logFC against logCPM with top FDR hairpins highlighted"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/plotSMEAR.R"       

rule camera:
    input:
        rds=["results/contrasts_matrix.rds", "results/model_matrix.rds", 
        "results/estimateDisp.rds"]
    output:
        tsv="results/{contrast}-camera.tsv",
        rds="results/{contrast}-camera.rds"
    params:
        contrast=get_contrast
    log:
        out = "logs/{contrast}-camera.out",
        err = "logs/{contrast}-camera.err"
    message:
        "Gene Set Test"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/camera.R"
 
rule stouffers:
    input:
        rds=["results/estimateDisp.rds", 
        "results/{contrast}-glmLRT.rds"]
    output:
        tsv="results/{contrast}-stouffers.tsv",
        rds="results/{contrast}-stouffers.rds"
    log:
        out = "logs/{contrast}-stouffers.out",
        err = "logs/{contrast}-stouffers.err"
    message:
        "Stouffer's method"
    conda:
        "../envs/stouffers.yaml"
    script:
        "../scripts/stouffers.R"

rule combinded_logFC:
    input:
        rds=["results/estimateDisp.rds", 
        "results/{contrast}-glmLRT.rds"]
    output:
        tsv="results/{contrast}-combinded_logFC.tsv",
        rds="results/{contrast}-combinded_logFC.rds"
    log:
        out = "logs/{contrast}-combinded_logFC.out",
        err = "logs/{contrast}-combinded_logFC.err"
    message:
        "Combinded logFC per gene"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/combinded_logFC.R"
  