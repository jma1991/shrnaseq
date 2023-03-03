# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

rule plotMDS:
    input:
        rds = "results/calcNormFactors.rds"
    output:
        png = "results/plotMDS.png"
    params:
        group = "Condition"
    log:
        out = "logs/plotMDS.png",
        err = "logs/plotMDS.png"
    message:
        "Plot Multi-Dimensional Scaling"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plotMDS.R"

rule plotPCA:
    input:
        rds = "results/calcNormFactors.rds"
    output:
        png = "results/plotPCA.png"
    params:
        group = "Condition"
    log:
        out = "logs/plotPCA.out",
        err = "logs/plotPCA.err"
    message:
        "Plot Principal Components Analysis"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plotPCA.R"

rule plotBCV:
    input:
        rds = "results/estimateDisp.rds"
    output:
        png = "results/plotBCV.png"
    log:
        out = "logs/plotBCV.out",
        err = "logs/plotBCV.err"
    message:
        "Plot Biological Coefficient of Variation"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plotBCV.R"

rule plotDist:
    input:
        rds = "results/calcNormFactors.rds"
    output:
        png = "results/plotDist.png"
    params:
        group = "Condition"
    log:
        out = "logs/plotDist.out",
        err = "logs/plotDist.err"   
    message:
        "Plot Sample Distance"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plotDist.R"
