analysis= function(input, output) {
    library(edgeR)
    x = processAmplicons(input$fastq, barcodefile=input$index, hairpinfile=input$hairpin, verbose=TRUE)
    saveRDS(x, file = output$rds)
}
analysis(snakemake@input, snakemake@output)
