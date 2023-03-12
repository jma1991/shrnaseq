# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule camera: 
    input:
        rds = ["results/estimateDisp.rds", "results/modelMatrix.rds", "results/makeContrasts.rds"]
    output:
        tsv = "results/{contrast}.camera.tsv"
    params:
        contrast = lambda wc: wc.contrast
    log:
        out = "logs/{contrast}.camera.out",
        err = "logs/{contrast}.camera.err"
    message:
        "Gene Set Test"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/camera.R"

rule combineTests: 
    input:
        tsv = "results/{contrast}.topTags.tsv"
    output:
        tsv = "results/{contrast}.combineTests.tsv"
    params:
        fdr = 0.05
    log:
        out = "logs/{contrast}.combineTests.out",
        err = "logs/{contrast}.combineTests.err"
    message:
        "Combine statistics across multiple tests"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/combineTests.R"
