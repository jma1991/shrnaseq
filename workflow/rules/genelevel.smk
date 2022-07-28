
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
        "../envs/edger.yaml"
    script:
        "../scripts/camera.R"
 
rule stouffers:
    input:
        rds=["results/estimateDisp.rds", 
        "results/{contrast}-glmLRT.rds"]
    output:
        tsv="results/{contrast}-stouffers.tsv",
        rds="results/{contrast}-stouffers.rds"
    log:
        out = "logs/{contrast}-stouffers.out",
        err = "logs/{contrast}-stouffers.err"
    message:
        "Stouffer's method"
    conda:
        "../envs/stouffers.yaml"
    script:
        "../scripts/stouffers.R"

rule combinded_logFC:
    input:
        rds=["results/estimateDisp.rds", 
        "results/{contrast}-glmLRT.rds"]
    output:
        tsv="results/{contrast}-combinded_logFC.tsv",
        rds="results/{contrast}-combinded_logFC.rds"
    log:
        out = "logs/{contrast}-combinded_logFC.out",
        err = "logs/{contrast}-combinded_logFC.err"
    message:
        "Combinded logFC per gene"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/combinded_logFC.R"
  