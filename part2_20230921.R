#20230921_study script part2

install.packages('dplyr')
install.packages('Seurat')

library(dplyr)
library(Seurat)
library(patchwork)
pbmc_small
pbmc <- pbmc_small

pbmc

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:17, cells = 500, balanced = TRUE)

ElbowPlot(pbmc)

pbmc
