# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

def get_final_output():
    output = [
        "results/processAmplicons.rds",
        "results/filterAmplicons.rds",
        "results/calcNormFactors.rds",
        "results/removeBatchEffect.rds",
        "results/modelMatrix.rds",
        "results/makeContrasts.rds",
        "results/glmFit.rds",
        "results/plotMDS.png",
        "results/plotPCA.png",
        "results/plotBCV.png",
        "results/plotDist.png"
        # "plots/counts-index-guideRNAs.png",
        # "plots/MDS-plot.png",
        # "plots/corrected-MDS-plot.png",
        # "plots/BCV-plot.png",
        # "plots/PCA-plot.png",
        # "plots/corrected-PCA-plot.png",
        # "plots/sample-dist-heatmap.png",
        # "plots/corrected-sample-dist-heatmap.png",
    ]

    contrasts = config["contrasts"]
    for contrast in contrasts:
        output.append(f"results/{contrast}.glmLRT.rds")
        output.append(f"results/{contrast}.topTags.tsv")
        output.append(f"results/{contrast}.camera.tsv")
        # output.append(f"plots/{contrast}-expression-heatmap.png")
        # output.append(f"plots/{contrast}-volcano-plot.png")
        # output.append(f"plots/{contrast}-guideRNA-histogram.png")
        # output.append(f"results/{contrast}-top-ranked-guideRNAs.tsv")
        # output.append(f"results/{contrast}-FDR-sig-guideRNAs.tsv")
        # output.append(f"results/{contrast}-FDR_guideRNAs.rds")
        # output.append(f"plots/{contrast}-plotSmear.png")
        # output.append(f"plots/{contrast}-camerarank.png")
        # output.append(f"plots/{contrast}-corrected-expression-heatmap.png")
        # output.append(f"results/{contrast}-gene-level.tsv")
        # output.append(f"results/{contrast}-gene-level.rds")
        # output.append(f"plots/{contrast}-generank.png")
        # output.append(f"results/{contrast}-goana.tsv")
        # output.append(f"results/{contrast}-goana.rds")
        # output.append(f"results/{contrast}-top_goana.tsv")
        # output.append(f"results/{contrast}-kegg.tsv")
        # output.append(f"results/{contrast}-kegg.rds")
        # output.append(f"results/{contrast}-top_kegg.tsv")
        # output.append(f"results/{contrast}-shinydata.rds")

    return output



def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]
