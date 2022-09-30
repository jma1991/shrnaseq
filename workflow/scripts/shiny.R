analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    x=readRDS(input$rds[1])
    corrected=readRDS(input$rds[2])
    des=readRDS(input$rds[3])
    matrix=readRDS(input$rds[4])
    xglm=readRDS(input$rds[5])
    lrt=readRDS(input$rds[6])
    assign(paste0(params$contrast, "_lrt"), lrt)

    go=readRDS(input$rds[7])
    assign(paste0(params$contrast, "_go"), go)
    kegg=readRDS(input$rds[8])
    assign(paste0(params$contrast, "_kegg"), kegg)
    camera=read.table(input$tsv[1], header=T)
    assign(paste0(params$contrast, "_camera"), camera)
    genelevel=read.table(input$tsv[2], header=T, sep=",")
    assign(paste0(params$contrast, "_genelevel"), genelevel)
    
    list=ls()
    list=list[! list %in% c("input", "output", "params", "log", "out", "err", "lrt", "go", "kegg", "camera", "genelevel")]

    list=mget(list)

    saveRDS(list, file = output$rds)
}


analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)