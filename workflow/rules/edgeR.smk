rule processAmplicons:
    input: 
        index=config["samples"],
        hairpins=config["hairpins"],
        fastq=config["fastq"]
    output: 
        rds="results/processAmplicons.rds"
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
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/filter_hairpins.R"

rule MDS_plot:
    input:
        rds="results/filter_hairpins.rds"
    output:
        plot="plots/MDS-plot.png"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/MDS_plot.R"   

rule diff_rep_analysis:
    input:
        rds="results/filter_hairpins.rds"
    output:
        Rdata="results/diff_rep_analysis.Rdata"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/diff_rep_analysis.R" 

rule BVC_plot:
    input:
        Rdata="results/diff_rep_analysis.Rdata"
    output:
        plot="plots/BCV-plot.png"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/BVC_plot.R" 
    
rule glmFit:
    input:
        Rdata="results/diff_rep_analysis.Rdata"
    output:
        rds="results/glmFit.rds"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmFit.R" 

rule glmLRT:
    input:
        rds="results/glmFit.rds"
    output:
        rds="results/glmLRT.rds"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmLRT.R" 

rule top_hairpins:
    input:
        rds="results/glmLRT.rds"
    output:
        txt="results/top-ranked-hairpins.txt"
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
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/FDR_hairpins.R"

rule plotSMEAR:
    input:
        rds=["results/glmLRT.rds", "results/FDR_hairpins.rds"]
    output:
        plot="plots/plotSmear.png"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/plotSMEAR.R"       


