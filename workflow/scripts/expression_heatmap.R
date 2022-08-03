hairpins_to_genes=function(mat,x){

  hairpinlist=list()
  unq = unique(x$genes$Gene)
  unq = unq[!is.na(unq)]
  
  for (i in unq) {
    sel = x$genes$Gene == i & !is.na(x$genes$Gene)
    hairpinlist[[i]] = x$genes$ID[c(which(sel))]
  }
  
  vector=vector()
  mat2=vector()
  for (i in 1:ncol(mat)) {
    for (i2 in hairpinlist) {
      vector=rbind(vector,mean(mat[c(i2),i]))
    }
    mat2=cbind(mat2, vector)
    vector=NULL
  }
  rownames(mat2)=unique(x$genes$Gene)
  colnames(mat2)=colnames(mat)
  return(mat2)
}

analysis=function(input, output, params, log) {
    
    #Log 
    out <- file(log$out, open = "wt")

    err <- file(log$err, open = "wt")

    sink(out, type = "output")

    sink(err, type = "message")

    #Script
    library(RColorBrewer)
    library(pheatmap)
    library(edgeR)
    
    lrt=readRDS(input$rds[1])
    top2 = topTags(lrt, n=Inf)

    x=readRDS(input$rds[2])
    mat <- cpm(x$counts, log=TRUE, prior.count = 1)
    colnames(mat) <- paste(x$samples$group, x$samples$Replicate, sep = " - " )
    #mat=hairpins_to_genes(mat,x)
    #selY <- top2$table$Gene[abs(top2$table$logFC)>1.5]
    selY <- rownames(top2$table)[abs(top2$table$logFC)>params$FC]
    
    mat <- subset(mat, rownames(mat) %in% selY)
    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

    png(output$plot[1], width=2000, height=1800, res=400)
    pheatmap(mat, col = colors, main="Differential expression across the groups (logCPM)")
    dev.off()

    #batch corrected
    mat=readRDS(input$rds[3])
    colnames(mat)= paste(x$samples$group, x$samples$Replicate, sep = " - " )
    #matrix=hairpins_to_genes(mat,x)
    selY <- rownames(top2$table)[abs(top2$table$logFC)>params$FC]
    mat = subset(mat, rownames(mat) %in% selY)
    colors <- colorRampPalette( brewer.pal(9, "Blues") )(255)

    png(output$plot[2], width=2000, height=1800, res=400)
    pheatmap(mat, col = colors, main="Batch corrected differential \n expression across the groups (logCPM)")
     dev.off()
}
  
analysis(snakemake@input, snakemake@output, snakemake@params, snakemake@log)