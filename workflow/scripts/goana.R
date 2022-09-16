analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(limma)
    library(params$organism, lib.loc=dirname(input$pkg), character.only = TRUE)

    matrix=readRDS(input$rds[1])
    lrt=readRDS(input$rds[2])

    org <- params$organism
    obj <- getFromNamespace(org, org)


    row.names(lrt) <- keys(obj)[1:nrow(lrt)]
    go <- goana(lrt, con=matrix[, params$contrast], FDR=params$threshold, 
    species = strsplit(params$organism, ".", fixed = TRUE)[[1]][2])
  
    write.table(go, file = output$tsv, quote = FALSE, sep = '\t', col.names = NA)
    saveRDS(go,file=output$rds)
}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)