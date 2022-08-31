analysis= function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    print(params$hairpinBeforeBarcode)
    x = processAmplicons(c(input$fastq), 
                        readfile2=params$readfile2,
                        barcodefile=input$index, 
                        hairpinfile=input$hairpin, 
                        verbose=TRUE,
                        hairpinBeforeBarcode=params$hairpinBeforeBarcode
                        )
    saveRDS(x, file = output$rds)
}
analysis(snakemake@input, snakemake@output,snakemake@params,snakemake@log)
