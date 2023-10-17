#part3_20231017
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/cloud/project/filtered_gene_bc_matrices/hg19/")
pbmc.data
# Initialize the Seurat object with the raw (non-normalized data).
#ppt
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_study", min.cells = 3, min.features = 200)
pbmc@meta.data %>% head
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_study")

pbmc@meta.data %>% dim

#MT-DNA percentage check
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#ppt next slide



# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2


#mitochondria contaminated cell filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#normalization 
pbmc <- NormalizeData(pbmc)


#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #top 2000 gene selection, variance stabilizing transformation

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(pbmc)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

#memory..

all.genes <- rownames(pbmc)
select <- all.genes[1:5000] 
pbmc <- ScaleData(pbmc, features = select)#

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

rm(pbmc.data)#memory save...


##############################3rd study start

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)

pbmc
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc #umap label+
DimPlot(pbmc, reduction = 'umap')



# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

cluster3vs1_5.markers <- FindMarkers(pbmc, ident.1 = 3, ident.2 = c(1, 5), min.pct = 0.25)
head(cluster3vs1_5.markers)

pbmc@meta.data %>% head


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% dim
pbmc.markers %>% filter(cluster == 0) %>% dim



#clustering 조정 
head(Idents(pbmc), 5) 

#pbmc <- FindNeighbors(pbmc, dims = 1:10)
#pbmc <- FindClusters(pbmc, resolution = 0.5)
#pbmc <- RunUMAP(pbmc, dims = 1:10)
#DimPlot(pbmc, reduction = 'umap')
head(pbmc@meta.data)
head(Idents(pbmc), 5)


#resolution 조정 
pbmc <- FindClusters(pbmc, resolution = 0.1)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')
pbmc <- FindClusters(pbmc, resolution = 1)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')

#dimension 조정
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:2)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:2)
head(pbmc@meta.data)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')




pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:20)
head(pbmc@meta.data)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')




pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:5)
head(pbmc@meta.data)
head(Idents(pbmc), 5)
DimPlot(pbmc, reduction = 'umap')






#re
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

head(Idents(pbmc), 5)


pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = 'umap')


VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
