# Authors: James Ashmore, Claire Prince
# Copyright: Copyright 2023, Zifo Technologies Ltd.
# Email: james.ashmore@zifornd.com
# License: MIT

def get_final_output(config):

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
    ]

    contrasts = config["contrasts"]
    for contrast in contrasts:
        output.append(f"results/{contrast}.glmLRT.rds")
        output.append(f"results/{contrast}.topTags.tsv")
        output.append(f"results/{contrast}.combineTests.tsv")
        output.append(f"results/{contrast}.camera.tsv")
        output.append(f"results/{contrast}.goana.tsv")
        output.append(f"results/{contrast}.kegga.tsv")

    return output

def get_contrast(wildcards):
    return config["contrasts"][wildcards.contrast]
