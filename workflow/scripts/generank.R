 analysis=function(input, output, params, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(edgeR)
    library(ggplot2)
    genelevel <- read.table(input$tsv, header=T, sep=",")

    colnames(genelevel) <- c("Gene", "nGuides", "Mean logFC", "IQR logFC", "Direction mean logFC",	"Direction smallest Pvalue", "Stouffer's Pvalue", "Stouffer's FDR")
    print(head(genelevel))
    genelevel$Rank <- rank(genelevel$"Mean logFC")
    print(head(genelevel))
    plt=ggplot(genelevel, aes(x=Rank, y=`Mean logFC`, label=Gene)) +
      geom_point(aes(size=nGuides), color="#b8dbcc") +
      #geom_text(data=genelevel[genelevel$`Mean logFC`>params$FC,], check_overlap = TRUE, nudge_x = max(genelevel$Rank)*0.05) +
      #geom_text(data=genelevel[genelevel$`Mean logFC`<(-(params$FC)),], check_overlap = TRUE, nudge_x = max(genelevel$Rank)*0.05) +
      labs(title= "Gene level rank plot") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_classic()
    png(output$plot, width=3000, height=3500, res=400)
    print(plt)
    dev.off()

}

analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)