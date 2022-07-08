analysis=function(input, output) {
    library(edgeR)
    load(input$Rdata)
    png(output$plot, width=2000, height=2000, res=400)
        plotBCV(xglm, main="Another small screen: BCV Plot")
    dev.off()
}

analysis(snakemake@input, snakemake@output)