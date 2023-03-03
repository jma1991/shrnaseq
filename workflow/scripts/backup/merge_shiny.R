analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    data=list()
    for (i in 1:length(input$rds)) {
        inputlist=readRDS(input$rds[i])
        data=c(data, inputlist)

    }
    saveRDS(data, file = output$rds)
}


analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)