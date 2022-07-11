analysis=function(input, output) {
    library(edgeR)
    xglm=readRDS(input$rds)
    png(output$plot, width=2000, height=2000, res=400)
        plotBCV(xglm, main="Another small screen: BCV Plot")
    dev.off()
}

analysis(snakemake@input, snakemake@output)