rule shiny:
    input:
        rds=["results/filter_guideRNAs.rds",
        "results/corrected_counts.rds",
        "results/model_matrix.rds",
        "results/contrasts_matrix.rds",
        "results/estimateDisp.rds"]
    output:
        rds="results/shinydata.rds"
    log:
        out = "logs/shiny.out",
        err = "logs/shiny.err"
    message:
        "Data for shiny app"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/shiny.R"

rule shiny_contrasts:
    input:
        rds=["results/{contrast}-glmLRT.rds",
        "results/{contrast}-goana.rds",
        "results/{contrast}-kegg.rds"],
        tsv=["results/{contrast}-camera.tsv",
        "results/{contrast}-gene-level.tsv"]
    output:
        rds="results/{contrast}-shinydata.rds"
    params:
        contrast=get_contrast
    log:
        out = "logs/{contrast}-shiny.out",
        err = "logs/{contrast}-shiny.err"
    message:
        "Contrast data for shiny app"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/shiny_contrasts.R"

rule merge_shiny:
    input:
        rds=['results/shinydata.rds',
        expand('results/{i}-shinydata.rds', i=config["contrast"])]
    output:
        rds="results/shiny.rds"
    log:
        out = "logs/merged-shiny.out",
        err = "logs/merged-shiny.err"
    message:
        "Merge shiny files"
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/merge_shiny.R"