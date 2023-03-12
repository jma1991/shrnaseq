# Author: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule goana:
    input:
        tsv = "results/{contrast}.combineTests.tsv",
        pkg = expand("resources/bioconductor/organism/lib/R/library/{package}", package = config["organism"])
    output:
        dir = directory("results/{contrast}.goana")
    params:
        FDR = 1.1,
        organism = config["organism"]
    log:
        out = "logs/{contrast}.goana.out",
        err = "logs/{contrast}.goana.err"
    message:
        "Test over-representation of GO terms for contrast: {wildcards.contrast}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/goana.R"

# rule topgo:
#     input:
#         tsv = "results/{contrast}.goana.tsv"
#     output:
#         pdf = "results/{contrast}.topgo.pdf"
#     params:
#         number = 20
#     log:
#         out = "logs/{contrast}.topgo.out",
#         err = "logs/{contrast}.topgo.err"
#     message:
#         "Plot top GO terms for contrast: {wildcards.contrast}"
#     conda:
#         "../envs/environment.yaml"
#     script:
#         "../scripts/topgo.R"

rule kegga:
    input:
        tsv = "results/{contrast}.combineTests.tsv",
        pkg = expand("resources/bioconductor/organism/lib/R/library/{package}", package = config["organism"])
    output:
        dir = directory("results/{contrast}.kegga")
    params:
        FDR = 1.1,
        organism = config["organism"]
    log:
        out = "logs/{contrast}.kegga.out",
        err = "logs/{contrast}.kegga.err"
    message:
        "Test over-representation of KEGG pathways for contrast: {wildcards.contrast}"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/kegga.R"

# rule topkegg:
#     input:
#         tsv = "results/{contrast}.kegga.tsv"
#     output:
#         pdf = "results/{contrast}.topkegg.pdf"
#     params:
#         number = 20
#     log:
#         out = "logs/{contrast}.topkegg.out",
#         err = "logs/{contrast}.topkegg.err"
#     message:
#         "Plot top KEGG pathways for contrast: {wildcards.contrast}"
#     conda:
#         "../envs/environment.yaml"
#     script:
#         "../scripts/topkegg.R"
