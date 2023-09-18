#20230921_study script part1

install.packages('dplyr')
install.packages('Seurat')

library(dplyr)
library(Seurat)
library(patchwork)

#file download
url <- 'https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
destfile <- '/cloud/project/pbmc3k_filtered_gene_bc_matrices.tar.gz'
download.file(url, destfile)

#terminal part
/cloud/project$ pwd
/cloud/project$ ls
pbmc3k_filtered_gene_bc_matrices.tar.gz  project.Rproj
/cloud/project$ tar -xf pbmc3k_filtered_gene_bc_matrices.tar.gz 
#################################################################

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/cloud/project/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
#ppt
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_study", min.cells = 3, min.features = 200)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc_study")

pbmc@meta.data %>% head

#MT-DNA percentage check
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#ppt next slide



# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#mitochondria contaminated cell filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#normalization 
pbmc <- NormalizeData(pbmc)


#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #top 2000 gene selection, variance stabilizing transformation

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#memory....

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
