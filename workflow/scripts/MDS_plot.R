mdsplt <- function(x, mat, mdstitle) {
        
    var <- matrixStats::rowVars(mat)
    num <- min(500, length(var))
    ind <- order(var, decreasing = TRUE)[seq_len(num)]
    dst <- dist(t(mat[ind, ]))
    mds <- cmdscale(as.matrix(dst))
    dat <- data.frame(
    MD1 = mds[, 1], 
    MD2 = mds[, 2], 
    group = x$samples$group
    )

    ggplot(dat, aes(MD1, MD2, colour = group)) + 
    geom_point(size = 3) + 
    labs(x = "MDS 1", y = "MDS 2", colour = "") + 
    labs(title=mdstitle) +  
    scale_color_brewer(palette = "Set3") +
    theme_classic()

}

analysis=function(input, output, log) {
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script 
    library(edgeR)
    library(RColorBrewer)
    library(ggplot2)
    x=readRDS(input$rds[1])
    mat <- cpm(x$counts, log = TRUE, prior.count = 1)
    png(output$plot[1], width=2500, height=1800, res=400)
    par(mar=c(5,5,3,7), xpd=TRUE)
    print(mdsplt(x, mat, "Multidimensional plot"))
    dev.off()
    
    ##batch corrected
    mat=readRDS(input$rds[2])

    png(output$plot[2], width=2500, height=1800, res=400)
    par(mar=c(5,5,3,7), xpd=TRUE)
    print(mdsplt(x, mat, "Batch corrected multidimensional plot"))
    dev.off()
}

analysis(snakemake@input, snakemake@output, snakemake@log)