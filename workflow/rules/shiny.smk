rule shiny:
    input:
        rds=["results/filter_hairpins.rds",
        "results/corrected_counts.rds",
        "results/model_matrix.rds",
        "results/contrasts_matrix.rds",
        "results/estimateDisp.rds",
        "results/{contrast}-glmLRT.rds",
        "results/{contrast}-FDR_hairpins.rds"]
    output:
        rdata="results/{contrast}-shiny.Rdata"
    log:
        out = "logs/{contrast}-shiny.out",
        err = "logs/{contrast}-shiny.err"
    message:
        "Data for shiny app"
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/shiny.R"
