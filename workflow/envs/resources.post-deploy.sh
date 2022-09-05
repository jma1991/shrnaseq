Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
Rscript -e 'BiocManager::install(version = "3.15")'
Rscript -e 'BiocManager::install(c("BiocGenerics", "S4Vectors", "IRanges", "AnnotationDbi"), lib="resources/bioconductor/organism/lib/R/library/")'
Rscript -e 'BiocManager::install("AnnotationDbi", lib="resources/bioconductor/organism/lib/R/library/")'
