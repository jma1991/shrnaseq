# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule modelMatrix:
    input:
        rds = "results/calcNormFactors.rds"
    output:
        rds = "results/modelMatrix.rds"
    log:
        out = "logs/modelMatrix.out",
        err = "logs/modelMatrix.err"
    message:
        "Construct Design Matrix"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/modelMatrix.R"

rule makeContrasts:
    input:
        rds = "results/calcNormFactors.rds"
    output:
        rds = "results/makeContrasts.rds"
    log:
        out = "logs/makeContrasts.out",
        err = "logs/makeContrasts.err"
    message:
        "Construct Contrasts Matrix"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/makeContrasts.R"

rule estimateDisp:
    input:
        rds = ["results/calcNormFactors.rds", "results/modelMatrix.rds"]
    output:
        rds = "results/estimateDisp.rds"
    log:
        out = "logs/estimateDisp.out",
        err = "logs/estimateDisp.err"
    message: 
        "Estimate Dispersions"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/estimateDisp.R" 

rule glmFit:
    input:
        rds = ["results/estimateDisp.rds", "results/modelMatrix.rds"]
    output:
        rds = "results/glmFit.rds"
    log:
        out = "logs/glmFit.out",
        err = "logs/glmFit.err" 
    message:
        "Fit NB-GLM"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/glmFit.R" 

rule glmLRT:
    input:
        rds = ["results/glmFit.rds", "results/makeContrasts.rds"]
    output:
        rds = "results/{contrast}.glmLRT.rds"
    params:
        contrast = lambda wildcards: wildcards.contrast
    log:
        out = "logs/{contrast}.glmLRT.out",
        err = "logs/{contrast}.glmLRT.err" 
    message:
        "Perform likelihood ratio test on GLM"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/glmLRT.R" 

rule topTags:
    input:
        rds = "results/{contrast}.glmLRT.rds"
    output:
        tsv = "results/{contrast}.topTags.tsv"
    log:
        out = "logs/{contrast}.topTags.out",
        err = "logs/{contrast}.topTags.err"
    message:
        "Table of the Top Differentially Expressed Tags"
    conda:
        "../envs/bioconductor-edger.yaml"
    script:
        "../scripts/topTags.R"
