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

# rule camerarank: 
#     input:
#         tsv="results/{contrast}-camera.tsv"
#     output:
#         plot="plots/{contrast}-camerarank.png"
#     params:
#         contrast=get_contrast,
#         FDR=config["FDR"]
#     log:
#         out = "logs/{contrast}-camerarank.out",
#         err = "logs/{contrast}-camerarank.err"
#     message:
#         "Camera rank plot"
#     conda:
#         "../envs/plots.yaml"
#     script:
#         "../scripts/camerarank.R"

# rule gene_level:
#     input:
#         rds="results/{contrast}-glmLRT.rds"
#     output:
#         tsv="results/{contrast}-gene-level.tsv",
#         rds="results/{contrast}-gene-level.rds"
#     log:
#         out = "logs/{contrast}-gene-level.out",
#         err = "logs/{contrast}-gene-level.err"
#     message:
#         "Gene level analysis"
#     conda:
#         "../envs/analysis.yaml"
#     script:
#         "../scripts/gene_level.R"

# rule generank: 
#     input:
#         tsv="results/{contrast}-gene-level.tsv"
#     output:
#         plot="plots/{contrast}-generank.png"
#     params:
#         contrast=get_contrast,
#         FC=config["FC"]
#     log:
#         out = "logs/{contrast}-generank.out",
#         err = "logs/{contrast}-generank.err"
#     message:
#         "Gene rank plot"
#     conda:
#         "../envs/plots.yaml"
#     script:
#         "../scripts/generank.R"
  