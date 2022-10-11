---
title: "single_cell"
output: html_document
---

Set the working directory.

```{r}
setwd("C:\\Users\\proca\\OneDrive\\Desktop\\vignette_scrna")
```

# Introduction

Attach the required libraries.

```{r warning=FALSE, message=FALSE}
library(Seurat)
library(dplyr)
library(patchwork)
```

Load the data.

```{r}
data <- get(load('SRA667466_SRS3060008.sparse.RData'))
```

# Pre-processing

Remove ENS part.

```{r}
rownames(data) <- gsub("\\_ENS.*", "", rownames(data))
head(rownames(data))
```
```{r}
## [1] "00R_AC107638.2" "0610009O20Rik"  "1110020A21Rik"  "1600014C23Rik" 
## [5] "1700007L15Rik"  "1700015F17Rik"
```

Create the Seurat Object.

```{r}
dorsal <- CreateSeuratObject(counts = data, project = "dorsal", min.cells = 3, min.features = 200)
dorsal
```
```{r}
## An object of class Seurat 
## 23577 features across 7546 samples within 1 assay 
## Active assay: RNA (23577 features, 0 variable features)
```

Calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of genes. Data are from mouse, we identify the mitochondrial genes starting with "mt-".

```{r}
dorsal[["percent.mt"]] <- PercentageFeatureSet(dorsal, pattern = "^mt-")
```

Show QC metrics for the first 5 cells.

```{r}
head(dorsal@meta.data, 5)
```
```{r}
##                orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACATACACGCTA     dorsal       1228          855   1.384365
## AAACATACGCGAAG     dorsal       7139         2545   1.470794
## AAACATACTGAGAA     dorsal       1004          732   1.593625
## AAACATACTGCTGA     dorsal        713          561   2.384292
## AAACATTGCATACG     dorsal       1129          781   2.302923
```
Violin plot.

```{r}
VlnPlot(dorsal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
```

```{r}
plot1 <- FeatureScatter(dorsal, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dorsal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

Final data. Looking at the previous plot we decided to filter for nFeature_RNA > 200 and nFeature_RNA < 4500 and percent.mt < 5

```{r}
dorsal <- subset(dorsal, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)
dorsal
```
```{r}
## An object of class Seurat 
## 23577 features across 7127 samples within 1 assay 
## Active assay: RNA (23577 features, 0 variable features)
```

Final violin plot. 

```{r}
plot1 <- FeatureScatter(dorsal, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dorsal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

## Normalization

Normalize the data. Seurat normalizes the gene counts in three steps. First divides gene counts by the total counts for each cell, then multiplies it by a scale factor (10,000 by default), and finally log-transforms the result.

```{r}
dorsal <- NormalizeData(dorsal, normalization.method = "LogNormalize", scale.factor = 10000)
```

We restrict the gene set to the “most variable” genes, hence those genes with the highest cell-to-cell variation in the dataset. We keep the top 2000 variable genes.

```{r}
dorsal <- FindVariableFeatures(dorsal, selection.method = "vst", nfeatures = 2000)
```

##Detection of highly variable genes

Identify the 10 most highly variable genes.

```{r}
top10 <- head(VariableFeatures(dorsal), 10)
top10
```
```{r}
 ## [1] "Ptgds"  "Acta2"  "Lyz2"   "Pf4"    "Crip1"  "Gfap"   "Vtn"    "Sst"   
 ## [9] "Dcn"    "Cd209f"
 ```

Plot variable features with and without the labels

```{r}
plot1 <- VariableFeaturePlot(dorsal)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

## Gene scaling and PCA

Scale the data. We apply a recommended step by Seurat for 10X data. All log-normalized counts are transformed so that they have mean 0 and variance 1 across all cells, regardless of the count values (high or low).

```{r}
all.genes <- rownames(dorsal)
dorsal <- ScaleData(dorsal, features = all.genes)
```

When data are scaled, we can regress out from the results unwanted sources of variation. We have computed already the “percent.mt" and removed cells with "percent.mt" < 5. We can also regress out the effect of cell cycle.

```{r}
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dorsal <- CellCycleScoring(dorsal, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
```

We use dimensionality reduction with PCA and plot the result. In our case it is not needed to regress out the effect of cell cycle, since we can notice from the plot that there are not clusters.

```{r message=FALSE }
dorsal <- RunPCA(dorsal, features = VariableFeatures(object = dorsal))
DimPlot(dorsal)
```

Now we want to decide how many principal components keep for the clustering. We use the elbow plot because of its intuitivness. It is a ranking of principal components based on the percentage of variance explained by each one. 

```{r}
ElbowPlot(dorsal, ndims = 20)
```

We observe an elbow around 12. Although the plot stabilizes around 25, I decided to keep 12 as number of dimensions first because the clusterings with increased dimensions are similar to the one with 12 and then because KNN clustering is affected by "curse of dimensionality" so it was better keep a not much high number of dimensions.

# Clustering

Seurat applies modularity optimization techniques to cluster cells, that iteratively group cells together. 
The “optimal” number of clusters is computed automatically. After many attempts I decided to use 12 dimensions and a resolution of 0.4, this was the best clustering. 

```{r}
dorsal <- FindNeighbors(dorsal, dims = 1:12)
dorsal <- FindClusters(dorsal, resolution = 0.4)
```
```{r}
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

## Number of nodes: 7127
## Number of edges: 231237

## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9506
## Number of communities: 19
## Elapsed time: 0 seconds
```

Run the UMAP.

```{r}
dorsal <- RunUMAP(dorsal, dims = 1:12)
```

Plot the UMAP.

```{r}
DimPlot(dorsal, reduction = "umap")
```

Find the number of cells in each cluster.

```{r}
table(dorsal$RNA_snn_res.0.4)
```
```{r}
##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
## 1574  958  720  517  461  446  439  431  231  230  214  176  163  159  118  101 
##   16   17   18 
##   79   64   46 
```

# Cell type markers

Find the marker genes.Now it’s the moment to find the genes that make us capable of assigning a cellular type to each cluster. In order to identify the specificity of a gene for a cell it is necessary to retrieve information in PanglaoDB or in scientific literature.

```{r}
dorsal.markers <- FindAllMarkers(dorsal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

Save the object and load as a table for the future

```{r}
write.table(dorsal.markers, "dorsal.markers.txt", sep ="\t", row.names = TRUE, col.names = TRUE)
markers_table <- read.table("dorsal.markers.txt", header = TRUE)
```

Find the top 5 markers for each cluster

```{r}
top5 <- markers_table %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5
```
```{r}
## # A tibble: 95 x 7
## # Groups:   cluster [19]
##    p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene         
##    <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>        
##  1     0       3.87 0.923 0.147         0       0 Slc6a11      
##  2     0       3.74 0.809 0.12          0       0 RP23-168E14.7
##  3     0       3.64 0.752 0.091         0       0 Htra1        
##  4     0       3.62 0.963 0.283         0       0 Aldoc        
##  5     0       3.50 0.656 0.046         0       0 Cldn10       
##  6     0       3.78 0.999 0.158         0       1 Mal          
##  7     0       3.50 1     0.62          0       1 Plp1         
##  8     0       3.35 0.987 0.134         0       1 Trf          
##  9     0       3.12 1     0.154         0       1 Mog          
## 10     0       3.05 1     0.948         0       1 Fth1         
## # ... with 85 more rows
```
Use violin plots to show expression probability distributions across clusters.

```{r}
VlnPlot(dorsal, features = c("Slc6a11", "Mal","Meg3"))
```

```{r}
VlnPlot(dorsal, features = c("Cldn5", "Gad2","Mag"))
```

```{r}
VlnPlot(dorsal, features = c("Ctss", "Pdgfra","Enpp6"))
```

```{r}
VlnPlot(dorsal, features = c("Gfap", "Gad1","Rgs5"))
```

```{r}
VlnPlot(dorsal, features = c("Pcsk1n", "Acta2","Ptgds"))
```

```{r}
VlnPlot(dorsal, features = c("Mbp", "Syp"))
```

```{r}
VlnPlot(dorsal, features = c("Pf4","RP23-100C7.3"))
```

Use feature plots to show feature expression on the UMAP

```{r}
FeaturePlot(dorsal, features = c("Slc6a11","Mal","Meg3","Cldn5","Gad2","Mag","Ctss", "Pdgfra","Enpp6"))
```

```{r}
FeaturePlot(dorsal, features = c("Gfap", "Gad1","Rgs5","Pcsk1n", "Acta2","Ptgds", "Mbp","Syp","Pf4"                                           ,"RP23-100C7.3"))
```

Assign a cell type to each cluster and plot the final UMAP

```{r}
new.cluster.ids <- c("Astrocytes1", "Oligodendrocytes1", "Neurons1", "Endothelial cells", "Neurons2", 
                     "Oligodendrocytes2", "Microglia", "OPC", "Oligodendrocytes3", "Astrocytes2",
                     "Neurons3","Pericytes", "Neurons4", "Smooth muscle cells",
                     "Fibroblasts", "Oligodendrocytes4", "Neurons5", "Macrophages","Ependymal cells")
names(new.cluster.ids) <- levels(dorsal)
dorsal <- RenameIdents(dorsal, new.cluster.ids)
DimPlot(dorsal, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```










