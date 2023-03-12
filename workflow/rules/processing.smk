# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule processAmplicons: 
    input:
        readfile = config["readfile"],
        barcodefile = config["barcodefile"],
        hairpinfile = config["hairpinfile"]
    output: 
        rds = "results/processAmplicons.rds"
    params:
        hairpinBeforeBarcode = config["hairpinBeforeBarcode"]
    log:
        out = "logs/processAmplicons.out",
        err = "logs/processAmplicons.err"
    message:
        "Process Amplicons Data"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/processAmplicons.R"

rule filterAmplicons:
    input:
        rds = "results/processAmplicons.rds"
    output:
        rds = "results/filterAmplicons.rds"
    log:
        out = "logs/filterAmplicons.out",
        err = "logs/filterAmplicons.err"
    message:
        "Filter Amplicons By Expression Level"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/filterAmplicons.R"

rule calcNormFactors:
    input:
        rds = "results/filterAmplicons.rds"
    output:
        rds = "results/calcNormFactors.rds"
    log:
        out = "logs/calcNormFactors.out",
        err = "logs/calcNormFactors.err"
    message:
        "Library Size Normalization"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/calcNormFactors.R"

rule removeBatchEffect: 
    input:
        rds = "results/calcNormFactors.rds"
    output:
        rds = "results/removeBatchEffect.rds"
    log:
        out = "logs/removeBatchEffect.out",
        err = "logs/removeBatchEffect.err" 
    message: 
        "Remove Batch Effect"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/removeBatchEffect.R" 
