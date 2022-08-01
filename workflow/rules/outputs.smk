def get_final_output():
    output = [
        "results/processAmplicons.rds",
        "plots/counts-index-hairpins.png",
        "results/corrected_counts.rds",
        "plots/MDS-plot.png",
        "plots/corrected-MDS-plot.png",
        "results/model_matrix.rds",
        "results/estimateDisp.rds",
        "plots/BCV-plot.png",
        "plots/PCA-plot.png",
        "plots/corrected-PCA-plot.png",
        "plots/sample-dist-heatmap.png",
        "plots/corrected-sample-dist-heatmap.png",
        "results/glmFit.rds"
    ]
    contrasts = config["contrast"]
    for contrast in contrasts:

        output.append(f"results/{contrast}-glmLRT.rds")
        output.append(f"plots/{contrast}-expression-heatmap.png")
        output.append(f"plots/{contrast}-volcano-plot.png")
        output.append(f"plots/{contrast}-hairpin-histogram.png")
        output.append(f"results/{contrast}-top-ranked-hairpins.tsv")
        output.append(f"results/{contrast}-FDR-sig-hairpins.tsv")
        output.append(f"results/{contrast}-FDR_hairpins.rds")
        output.append(f"plots/{contrast}-plotSmear.png")
        output.append(f"results/{contrast}-camera.tsv")
        output.append(f"results/{contrast}-camera.rds")
        output.append(f"plots/{contrast}-corrected-expression-heatmap.png")
        output.append(f"results/{contrast}-stouffers.tsv")
        output.append(f"results/{contrast}-stouffers.rds")
        output.append(f"results/{contrast}-combinded_logFC.tsv")
        output.append(f"results/{contrast}-combinded_logFC.rds")

    return output