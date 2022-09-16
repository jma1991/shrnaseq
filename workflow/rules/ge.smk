rule goana:
    input:
        rds = ["results/contrasts_matrix.rds", "results/{contrast}-glmLRT.rds"],
        pkg = expand("resources/bioconductor/organism/lib/R/library/{organism}", organism = config["organism"])
    output:
        tsv="results/{contrast}-goana.tsv",
        rds="results/{contrast}-goana.rds"
    params:
        contrast = get_contrast,
        threshold= config["FDR"],
        organism = config["organism"]
    log:
        out = "logs/{contrast}.goana.out",
        err = "logs/{contrast}.goana.err"
    message:
        "Test over-representation of GO terms for contrast"
    conda:
        "../envs/ge.yaml"
    script:
        "../scripts/goana.R"

rule top_goana:
    input:
        rds = "results/{contrast}-goana.rds"
    output:
        tsv="results/{contrast}-top_goana.tsv"
    log:
        out = "logs/{contrast}.top_goana.out",
        err = "logs/{contrast}.top_goana.err"
    message:
        "Test over-representation of GO terms for contrast"
    conda:
        "../envs/ge.yaml"
    script:
        "../scripts/top_goana.R"

rule kegg:
    input:
        rds = ["results/contrasts_matrix.rds", "results/{contrast}-glmLRT.rds"],
        pkg = expand("resources/bioconductor/organism/lib/R/library/{organism}", organism = config["organism"])
    output:
        tsv="results/{contrast}-kegg.tsv",
        rds="results/{contrast}-kegg.rds"
    params:
        contrast = get_contrast,
        threshold= config["FDR"],
        organism = config["organism"]
    log:
        out = "logs/{contrast}.kegg.out",
        err = "logs/{contrast}.kegg.err"
    message:
        "KEGG pathways"
    conda:
        "../envs/ge.yaml"
    script:
        "../scripts/kegg.R" 

rule top_kegg:
    input:
        rds ="results/{contrast}-kegg.rds"
    output:
        tsv="results/{contrast}-top_kegg.tsv"
    log:
        out = "logs/{contrast}.top_kegg.out",
        err = "logs/{contrast}.top_kegg.err"
    message:
        "KEGG pathways"
    conda:
        "../envs/ge.yaml"
    script:
        "../scripts/top_kegg.R" 