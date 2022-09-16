rule corrected_counts:
    input:
        rds="results/estimateDisp.rds"
    output:
        rds="results/corrected_counts.rds"
    log:
        out = "logs/corrected_counts.out",
        err = "logs/corrected_counts.err" 
    message: 
        "Remove batch effect"
    conda:
        "../envs/limma.yaml"
    script:
        "../scripts/corrected_counts.R" 

rule model_matrix:
    input:
        rds="results/norm.rds"
    output:
        rds="results/model_matrix.rds"
    log:
        out = "logs/model_matrix.out",
        err = "logs/model_matrix.err"
    message:
        "Generate model matrix for edgeR analysis"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/model_matrix.R"

rule contrasts_matrix:
    input:
        rds="results/model_matrix.rds"
    output:
        rds="results/contrasts_matrix.rds"
    log:
        out = "logs/contrasts_matrix.out",
        err = "logs/contrasts_matrix.err"
    message:
        "Generate matrix for group contrasts"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/contrasts_matrix.R"

rule estimateDisp:
    input:
        rds=["results/norm.rds", "results/model_matrix.rds"]
    output:
        rds="results/estimateDisp.rds"
    log:
        out = "logs/estimateDisp.out",
        err = "logs/estimateDisp.err"
    message: 
        "Perform differential representation analysis"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/estimateDisp.R" 

rule glmFit:
    input:
        rds=["results/model_matrix.rds", "results/estimateDisp.rds"]
    output:
        rds="results/glmFit.rds"
    log:
        out = "logs/glmFit.out",
        err = "logs/glmFit.err" 
    message:
        "Fit negative bionomial GLM"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmFit.R" 

rule glmLRT:
    input:
        rds=["results/contrasts_matrix.rds", "results/glmFit.rds"]
    output:
        rds="results/{contrast}-glmLRT.rds"
    params:
        contrast=get_contrast
    log:
        out = "logs/{contrast}-glmLRT.out",
        err = "logs/{contrast}-glmLRT.err" 
    message:
        "Perform likelihood ratio test on GLM"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/glmLRT.R" 

rule top_hairpins:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        tsv="results/{contrast}-top-ranked-hairpins.tsv"
    log:
        out = "logs/{contrast}-top-ranked-hairpins.out",
        err = "logs/{contrast}-top-ranked-hairpins.err"
    message:
        "Generate table of the top differentially expressed hairpins"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/top_hairpins.R"

rule FDR_hairpins:
    input:
        rds="results/{contrast}-glmLRT.rds"
    output:
        tsv="results/{contrast}-FDR-sig-hairpins.tsv",
        rds="results/{contrast}-FDR_hairpins.rds"
    params:
        threshold=config["FDR"]
    log:
        out = "logs/{contrast}-FDR-sig-hairpins.out",
        err = "logs/{contrast}-FDR-sig-hairpins.err"
    message:
        "Highlight and generate table of hairpins with FDR < 0.05"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/FDR_hairpins.R"