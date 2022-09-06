analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(params$organism, lib.loc="resources/bioconductor/organism/lib/R/library/", character.only = TRUE)

    matrix=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])
    row.names(lrt) <- keys(org.Mm.eg.db)[1:nrow(lrt)]
    go <- goana(lrt, con=matrix[, params$contrast], FDR=params$threshold,
    species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2])
    topgo <- topGO(go, sort="up")

    write.table(go, file = output$tsv, quote = FALSE, sep = '\t', col.names = NA)
}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)