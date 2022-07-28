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
