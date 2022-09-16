analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(limma)
    library(AnnotationDbi)
    library(params$organism, lib.loc=dirname(input$pkg), character.only = TRUE)

    matrix=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])
    org <- params$organism
    obj <- getFromNamespace(org, org)
    row.names(lrt) <- keys(obj)[1:nrow(lrt)]
    keg = kegga(lrt, species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2], 
    coef=matrix[, params$contrast],  FDR=params$threshold)
 
    topkegg <- topKEGG(keg, sort="up")

    write.table(keg, file = output$tsv, quote = FALSE, sep = '\t', col.names = NA)
    saveRDS(keg,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)