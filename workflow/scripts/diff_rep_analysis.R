analysis=function(input, output) {
    library(edgeR)
    x=readRDS(input$rds)
    des = model.matrix(~x$samples$group)
    xglm = estimateDisp(x, des)
    com_disp=sqrt(xglm$common.disp)
    print(com_disp)
    save(des,xglm,file=output$Rdata)
}

analysis(snakemake@input, snakemake@output)