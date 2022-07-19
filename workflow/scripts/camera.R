analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    matrix=readRDS(input$rds[1])
    des=readRDS(input$rds[2])
    xglm=readRDS(input$rds[3])

    genesymbols = xglm$genes$Gene
    genesymbollist = list()
    unq = unique(genesymbols)
    unq = unq[!is.na(unq)]

    for (i in unq) {
    sel = genesymbols == i & !is.na(genesymbols)
    if (sum(sel) > 3) 
        genesymbollist[[i]] =which(sel)
    }
    camera.res = camera(xglm, index = genesymbollist, des, contrast=matrix[, params$contrast])
    camera.res=camera.res[!is.na(camera.res$FDR),]
    camera.res$gene=rownames(camera.res)
    write.table(camera.res, output$tsv, quote=F, row.names=F)
    saveRDS(camera.res,file=output$rds)

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)