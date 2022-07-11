analysis=function(input, output) {
    library(edgeR)
    x=readRDS(input$rds[1])
    des=readRDS(input$rds[2])
    xglm = estimateDisp(x, des)
    com_disp=sqrt(xglm$common.disp)
    print(com_disp)
    saveRDS(xglm,file=output$rds)
}

analysis(snakemake@input, snakemake@output)