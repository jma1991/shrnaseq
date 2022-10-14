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
    
    data=ls()
    data=data[! data %in% c("input", "output", "params", "log", "out", "err")]

    data=mget(data)

    saveRDS(data, file = output$rds)
}


analysis(snakemake@input, snakemake@output, snakemake@log)