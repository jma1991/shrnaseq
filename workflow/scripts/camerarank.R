analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(ggplot2)
    camera <- read.table(input$tsv, header=T)
    colnames(camera) <- c("nGuides", "Direction", "Pvalue", "FDR", "Gene")
    camera$Rank <- rank(-log(camera$Pvalue))
    plt=ggplot(camera, aes(x=Rank, y=(-log(Pvalue)), label=Gene)) + 
    geom_point(aes(size=nGuides), color="#b8dbcc") +
    geom_text(data=camera[camera$FDR<params$FDR,], check_overlap = TRUE, nudge_x = max(camera$Rank)*0.05) +
    labs(title= "Camera P value rank plot") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_classic()
    png(output$plot, width=3000, height=3500, res=400)
    print(plt)
    dev.off()

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)