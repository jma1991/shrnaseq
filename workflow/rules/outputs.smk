def get_final_output():
    output = [
        "results/processAmplicons.rds",
        "plots/counts-index-hairpins.png",
        "plots/MDS-plot.png",
        "results/model_matrix.rds",
        "results/diff_rep_analysis.rds",
        "plots/BCV-plot.png",
        "plots/PCA-plot.png",
        "plots/sample-dist-heatmap.png",
        "results/glmFit.rds"
    ]
    contrasts = config["contrast"]
    for x in contrasts:

        output.append(f"results/{x}-glmLRT.rds")
        output.append(f"plots/{x}-expression-heatmap.png")
        output.append(f"plots/{x}-volcano-plot.png")
        output.append(f"plots/{x}-hairpin-histogram.png")
        output.append(f"results/{x}-top-ranked-hairpins.tsv")
        output.append(f"results/{x}-FDR-sig-hairpins.tsv")
        output.append(f"results/{x}-FDR_hairpins.rds")
        output.append(f"plots/{x}-plotSmear.png")

    return output