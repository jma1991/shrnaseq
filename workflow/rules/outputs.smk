def get_final_output():
    output = [
        directory(expand("resources/bioconductor/organism/lib/R/library/{organism}", organism = config["organism"])),
        "results/processAmplicons.rds",
        "plots/counts-index-hairpins.png",
        "results/filter_hairpins.rds",
        "results/norm.rds",
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
        output.append(f"results/{contrast}-gene-level.tsv")
        output.append(f"results/{contrast}-gene-level.rds")
        output.append(f"results/{contrast}-goana.tsv")
        output.append(f"results/{contrast}-goana.rds")
        output.append(f"results/{contrast}-top_goana.tsv")
        output.append(f"results/{contrast}-kegg.tsv")
        output.append(f"results/{contrast}-kegg.rds")
        output.append(f"results/{contrast}-top_kegg.tsv")
        output.append(f"results/{contrast}-shiny.RData")

    return output