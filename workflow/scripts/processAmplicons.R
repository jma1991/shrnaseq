analysis= function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    x = processAmplicons(input$fastq, barcodefile=input$index, hairpinfile=input$hairpin, verbose=TRUE)
    saveRDS(x, file = output$rds)
}
analysis(snakemake@input, snakemake@output, snakemake@log)
