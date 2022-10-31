
rule camera: 
    input:
        rds=["results/contrasts_matrix.rds", "results/model_matrix.rds", 
        "results/estimateDisp.rds"]
    output:
        tsv="results/{contrast}-camera.tsv",
        rds="results/{contrast}-camera.rds"
    params:
        contrast=get_contrast
    log:
        out = "logs/{contrast}-camera.out",
        err = "logs/{contrast}-camera.err"
    message:
        "Gene Set Test"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/camera.R"
 
rule gene_level:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        tsv="results/{contrast}-gene-level.tsv",
        rds="results/{contrast}-gene-level.rds"
    log:
        out = "logs/{contrast}-gene-level.out",
        err = "logs/{contrast}-gene-level.err"
    message:
        "Gene level analysis"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/gene_level.R"


  