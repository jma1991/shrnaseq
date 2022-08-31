rule goana:
    input:
        rds = ["results/contrasts_matrix.rds", "results/{contrast}-glmLRT.rds"]
    output:
        tsv="results/{contrast}-goana.tsv"
    params:
        contrast = get_contrast,
        threshold=config["FDR"],
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