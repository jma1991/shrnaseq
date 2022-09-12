analysis=function(input, output, log) {
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
    top2ids=readRDS(input$rds[7])
    go=readRDS(input$rds[8])
    kegg=readRDS(input$rds[9])

    save(list = ls(all.names = TRUE),file = output$rdata)
}


analysis(snakemake@input, snakemake@output, snakemake@log)