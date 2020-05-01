ssh nscc04-ib0

Andy human wk17 fetal kidney
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3073089
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3073088

/home/users/ntu/l160005/scratch/andy_kidneywk17
cds <- load_cellranger_data("/home/users/ntu/l160005/scratch/andy_kidneywk17")

# Download the files to desktop
# Rename to "matrix.mtz.gz", "features.tsv.gz" and "barcodes.tsv.gz"

R
library(dplyr)
library(Seurat)
library(patchwork)

list.files("/home/users/ntu/l160005/scratch/andy_kidneywk17/outs/filtered_gene_bc_matrices/GRCh37")
kidney.data1 <- Read10X(data.dir = "/home/users/ntu/l160005/scratch/andy_kidneywk17_1/outs/filtered_gene_bc_matrices/GRCh37")
kidney.data2 <- Read10X(data.dir = "/home/users/ntu/l160005/scratch/andy_kidneywk17_2/outs/filtered_gene_bc_matrices/GRCh37")
ls()

kidney1 <- CreateSeuratObject(counts = kidney.data1, project = "Kidney_Wk17_1")
kidney2 <- CreateSeuratObject(counts = kidney.data2, project = "Kidney_Wk17_2")

kidney.combined <- merge(kidney1, y = kidney2, add.cell.ids = c("Kidney_Wk17_1", "Kidney_Wk17_2"), project = "Kidney_Wk17")
kidney.combined

# Check if the combination was successful
head(colnames(kidney.combined))
table(kidney.combined$orig.ident)
# Kidney_Wk17_1 Kidney_Wk17_2
#         4166          3699

kidney.combined[["percent.mt"]] <- PercentageFeatureSet(kidney.combined, pattern = "^MT-")

tiff("qc_5.tiff", width=1200, height= 1000)
VlnPlot(kidney.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

kidney.filtered <- subset(kidney.combined, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 5 & nCount_RNA < 10000)

tiff("qc_passed.tiff", width=1200, height= 1000)
VlnPlot(kidney.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

saveRDS(kidney.filtered, "kidney_filtered.rds")

tiff("fscat_befQC.tiff", width=1200, height= 600)
plot1 <- FeatureScatter(kidney.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kidney.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

tiff("fscat_aftQC.tiff", width=1200, height= 600)
plot1 <- FeatureScatter(kidney.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kidney.filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# Normalize the data
kidney <- NormalizeData(kidney.filtered)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top <- head(VariableFeatures(kidney), 10)
# plot variable features with and without labels
tiff("top10", width=1600, height= 600)
plot1 <- VariableFeaturePlot(kidney)
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2
dev.off()

# Identify the 10 most highly variable genes
top <- head(VariableFeatures(kidney), 20)
# plot variable features with and without labels
tiff("top20", width=1600, height= 600)
plot1 <- VariableFeaturePlot(kidney)
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2
dev.off()

# Identify the 10 most highly variable genes
top <- head(VariableFeatures(kidney), 30)
# plot variable features with and without labels
tiff("top30", width=1600, height= 600)
plot1 <- VariableFeaturePlot(kidney)
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2
dev.off()

# Reload due to some technical issues and repeat processing
kidney <- readRDS("kidney_filtered.rds")
kidney <- NormalizeData(kidney)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all.genes <- rownames(kidney)
kidney <- ScaleData(kidney, features = all.genes)

# Perform linear dimensional reduction
kidney <- RunPCA(kidney, features = VariableFeatures(object = kidney))
# Examine and visualize PCA results a few different ways
print(kidney[["pca"]], dims = 1:5, nfeatures = 5)

tiff("vizdim.tiff", width=800, height= 600)
VizDimLoadings(kidney, dims = 1:2, reduction = "pca")
dev.off()

tiff("dimplot.tiff", width=800, height= 600)
DimPlot(kidney, reduction = "pca")
dev.off()

tiff("dimheatmap2.tiff", width=1200, height= 1800)
DimHeatmap(kidney, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

kidney <- JackStraw(kidney, num.replicate = 100)

saveRDS(kidney, "kidney_beforejackstraw.rds")
kidney <- readRDS("kidney_beforejackstraw.rds")

kidney <- ScoreJackStraw(kidney, dims = 1:20)

tiff("jackstraw.tiff", width=800, height= 600)
JackStrawPlot(kidney, dims = 1:20)
dev.off()

tiff("elbow.tiff", width=800, height= 600)
ElbowPlot(kidney)
dev.off()

# Find clusters
kidney <- FindNeighbors(kidney, dims = 1:20)
kidney <- FindClusters(kidney, resolution = 0.5)
head(Idents(kidney), 5)

# UMAP analysis
kidney <- RunUMAP(kidney, dims = 1:20)

# Plot UMAP 
tiff("umap.tiff", width=1200, height= 800)
DimPlot(kidney, reduction = "umap", label = TRUE)
dev.off()

saveRDS(kidney, file = "kidney_umapready.rds")

cluster1.markers <- FindMarkers(kidney, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
kidney.markers <- FindAllMarkers(kidney, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
kidney.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

       p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene
       <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
 1 0.             1.26 0.933 0.306 0.        0       NNAT
 2 0.             1.14 0.782 0.225 0.        1       CYP1B1
 3 0.             1.83 0.637 0.088 0.        2       SULT1E1
 4 1.80e-265      1.19 0.826 0.25  5.88e-261 3       LYPD1
 5 8.23e-198      1.61 0.945 0.388 2.69e-193 4       HMGB2
 6 0.             2.65 0.939 0.112 0.        5       TM4SF1
 7 3.22e- 78      2.30 0.83  0.451 1.06e- 73 6       HIST1H4C
 8 0.             3.46 0.875 0.035 0.        7       LUM
 9 8.65e-134      4.83 0.345 0.032 2.83e-129 8       REN
10 0.             1.96 0.659 0.03  0.        9       CLDN4
11 1.52e-233      2.56 0.935 0.068 4.99e-229 10      MAFB
12 0.             3.59 0.933 0.01  0.        11      SRGN
13 1.18e-206      2.87 0.947 0.066 3.87e-202 12      SERPINI1

cluster1.markers <- FindMarkers(kidney, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers)
ITM2C  0.905 1.0863933 0.810 0.994 0.838
NNAT   0.879 1.2610544 0.758 0.933 0.306
BST2   0.878 1.2497167 0.756 0.944 0.397
RPL3   0.877 0.4264968 0.754 1.000 1.000
RPS4X  0.867 0.4395273 0.734 0.999 1.000
CRABP2 0.830 1.0669068 0.660 0.892 0.410

tiff("dp_rep2.tiff", width=800, height= 600)
DotPlot(kidney, features = c("NNAT", "CYP1B1", "SULT1E1", "LYPD1", "HMGB2", "TM4SF1", "HIST1H4C", "LUM", "REN", "CLDN4", "MAFB", "SRGN", "SERPINI1")) + RotatedAxis()
dev.off()

tiff("featureplot2.tiff", width=1800, height= 1200)
FeaturePlot(kidney, features = c("NNAT", "CYP1B1", "SULT1E1", "LYPD1", "HMGB2", "TM4SF1", "HIST1H4C", "LUM", "REN", "CLDN4", "MAFB", "SRGN", "SERPINI1"))
dev.off()

tiff("top10_allclust.tiff", width=1800, height= 1200)
top10 <- kidney.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(kidney, features = top10$gene) + NoLegend()
dev.off()


# CLUSTER 5 - ENDOTHELIAL CELLS
tiff("vlnplot_clust5.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CD34","CLDN5","PLVAP","ENG"))
dev.off()

tiff("fplot_clust5.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("CD34","CLDN5","PLVAP","ENG"))
dev.off()

tiff("fplot_clust5.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("CD34","CLDN5","PLVAP","ENG","CD74","ECSCR","EGFL7","ESAM","SDPR"))
dev.off()

tiff("vlnplot_clust5.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CD34","CLDN5","PLVAP","ENG","CD74","ECSCR","EGFL7","ESAM","SDPR"))
dev.off()

tiff("fplot_clust5b.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("CD34","KDR","PLVAP","ENG","ESAM","ECSCR","EGFL7","SOX17","SDPR"))
dev.off()

tiff("vlnplot_clust5b.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CD34","KDR","PLVAP","ENG","ESAM","ECSCR","EGFL7","SOX17","SDPR"))
dev.off()

# CLUSTER 10 - PODOCYTES
cluster10.markers <- FindMarkers(kidney, ident.1 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster10.markers, n=40)

tiff("fplot_clust10.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("MAFB","NPHS2","PODXL","WT1","BCAM","CITED2","SYNPO","OLFM3","ASS1"))
dev.off()

tiff("vlnplot_clust10.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MAFB","NPHS2","PODXL","WT1","BCAM","CITED2","SYNPO","OLFM3","ASS1"))
dev.off()

tiff("fplot_clust10b.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("MAFB","NPHS2","PODXL","WT1","VEGFA","OLFM3","SYNPO","OLFM3","ASS1"))
dev.off()

tiff("vlnplot_clust10b.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MAFB","NPHS2","PODXL","WT1","VEGFA","OLFM3","SYNPO","OLFM3","ASS1"))
dev.off()

# SCREENING
tiff("vlnplot_clustx1.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("JAG1","PAX8","EMX2","LHX1","SYNPO","ALDH1A1","SLC12A1","SLC3A1","KRT8"))
dev.off()

tiff("vlnplot_clustx2.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CLDN5","ITGA8","SIX1","SIX2","SALL1","CITED1","BMP4","CDH1","CDH6"))
dev.off()

tiff("vlnplot_clustx3.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("SOX9","GATA3","LRP2","CUBN","SLC34A1","CLDN6","TCF21","OLFM3","ASS1"))
dev.off()


kidney <- readRDS("kidney_umapready.rds")
cluster0.markers <- FindMarkers(kidney, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n=30)

tiff("vlnplot_clustx4.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("ITM2C","NNAT","RPL3","EPCAM","PRDX2","RPS4X","CRABP2","DAPL1","PRDX2"))
dev.off()

# CLUSTER 0, 3 and 4 are NPCs
tiff("vlnplot_NPC.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("SIX1","SIX2","SALL1","CITED1","ITGA8","PAX8","EPCAM","NNAT","DAPL1"))
dev.off()

tiff("featureplot_NPC.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("SIX1","SIX2","SALL1","CITED1","ITGA8","PAX8","EPCAM","NNAT","DAPL1"))
dev.off()


tiff("vlnplot_clustx5pt.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("IGFBP7","TSPAN12","CUBN","SLC3A1","GLYATL1","CDH6","FXYD2","KRT18","CLU"))
dev.off()

tiff("vlnplot_clustx5dt.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("LIMCH1","MECOM","MAL","FXYD2","SLC12A1","ALDH1A1","ZNF44","CDH1","GATA3"))
dev.off()

tiff("featureplot_clustx5pt.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("IGFBP7","TSPAN12","CUBN","SLC3A1","GLYATL1","CDH6","FXYD2","KRT18","CLU"))
dev.off()

tiff("featureplot_clustx5dt.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("LIMCH1","MECOM","MAL","FXYD2","SLC12A1","ALDH1A1","ZNF44","CDH1","GATA3"))
dev.off()

cluster6.markers <- FindMarkers(kidney, ident.1 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster6.markers, n=40)

tiff("vlnplot_clust12.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("SERPINI1","IGF1","RSPO3","ISLR","DCN","CEBPD","PEAR1","IRF1","NR4A1"))
dev.off()

tiff("vlnplot_clust12b.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("TUBB3","MAP2","SOX2","GFAP","STMN1","HES1"))
dev.off()

tiff("vlnplot_clust12c.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CALB1","DCX","RBFOX3","S100B","PHOX2B"))
dev.off()


tiff("vlnplot_clust6a.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("TUBA1B","TUBB","KIAA0101","TYMS","HMGB2","NUSAP1","ZWINT","HMGB1","SMC2"))
dev.off()

tiff("vlnplot_clust6b.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CDK1","SMC4","CKS1B","HMGN2","PTN","CKS2","HIST1H4C","RARRES2","TOP2A"))
dev.off()

tiff("vlnplot_clust6c.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CLSPN","RAD51AP1","RANBP1","PRC1","CENPF","MAD2L1","NUCKS1","UBE2C","UBE2T"))
dev.off()

tiff("vlnplot_clust6d.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("PCNA","NUSAP1","TOP2A","ID1","ID2","ID3","PAX2","MEOX1","FOXD1"))
dev.off()

tiff("vlnplot_clust6e.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("AQP1","AQP2","SLC4A1","AQP6","SLC26A4","HMX2","SYT7","PARM1","ATP6V1B1"))
dev.off()

tiff("vlnplot_clust6f.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("SLC27A2","LRP2","UMOD","VIM","NRP1","S100A8","S100A9","PLAC8","S100A4"))
dev.off()

cluster11.markers <- FindMarkers(kidney, ident.1 = 11, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster11.markers, n=60)

tiff("vlnplot_clust11a.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("CD74","CXCR4","NEAT1","SAT1","CD53","GPX1","HCST","CD83"))
dev.off()

tiff("vlnplot_clust11b.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MEIS1","MEIS2","MEIS3","P4HB","PDGFRA","PDGFRB"))
dev.off()

tiff("vlnplot_clust11c.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("","","","","","","","",""))
dev.off()

tiff("vlnplot_clust9ub.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("ALDH1A1","CDH1","GATA3","SOX9","KRT8","FXYD2","JAG1","CALB1","AQP2"))
dev.off()

tiff("featureplot_clust9ub.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("ALDH1A1","CDH1","GATA3","SOX9","KRT8","FXYD2","JAG1","CALB1","AQP2"))
dev.off()

tiff("vlnplot_pt.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("NOTCH1","NOTCH2","HES1","NRARP","CTTNB"))
dev.off()

tiff("vlnplot_fib.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MEIS1", "MEIS2", "PDGFRA", "PDGFRB", "CD105", "CD29", "CD44", "CD90", "HSP47", "S100A4", "VIM", "SMA", "FOXD1", "FSP1"))
dev.off()

tiff("vlnplot_prolif.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("ZWINT", "PCNA", "NUSAP1", "HMGB2", "TYMS", "SMC2", "CDK1", "SMC4","CKS1B"))
dev.off()

tiff("vlnplot_prolif2.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("HIST1H4C", "TOP2A", "UBE2C", "MAD2L1", "CENPF", "RAD51AP1", "TUBA1B", "TUBB", "RANBP1"))
dev.off()

tiff("featureplot_prolif.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("ZWINT", "PCNA", "NUSAP1", "HMGB2", "TYMS", "SMC2", "CDK1", "SMC4","CKS1B"))
dev.off()

tiff("featureplot_prolif2.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("HIST1H4C", "TOP2A", "UBE2C", "MAD2L1", "CENPF", "RAD51AP1", "TUBA1B", "TUBB", "RANBP1"))
dev.off()

tiff("vlnplot_strom.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MEIS1", "MEIS2", "FOXD1", "PDGFRA", "PDGFRB", "WT1", "ISLR", "RSPO1", "RSPO2", "RSPO3", "IRF1", "PEAR1", "MEOX1", "MALAT1", "CTGF", "BMP2", "BMP4", "BMP7", "MORC1", "EMX2", "HNF1B", "HNF4A"))
dev.off()

tiff("vlnplot_strom.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("MEIS1", "MEIS2", "FOXD1", "PDGFRA", "PDGFRB", "WT1", "ISLR", "IRF1", "PEAR1"))
dev.off()

tiff("featureplot_strom.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("MEIS1", "MEIS2", "FOXD1", "PDGFRA", "PDGFRB", "WT1", "ISLR", "IRF1", "PEAR1"))
dev.off()

R
library(dplyr)
library(Seurat)
library(patchwork)
kidney <- readRDS("kidney_umapready.rds")

0 NPC 
1 Stromal I
2 Stromal II
3 Proximal tubule CDH1 PAX8
4 Proliferative NPCs
5 Endothelial
6 Proliferative stromal
7 Stromal III
8 Stromal IV
9 Collecting duct/Distal connecting tubule
10 Podocytes
11 Unidentified
12 Stromal V

tiff("vlnplot_clust11.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("HES1", "ID1", "CD74", "NRARP", "CXCR4", "CD53", "HCST", "CD83", "CD44"))
dev.off()

tiff("featureplot_clust11.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("HES1", "ID1", "CD74", "NRARP", "CXCR4", "CD53", "HCST", "CD83", "CD44"))
dev.off()

#### Code to rename cluster
new.cluster.ids <- c("NPC", "Stromal I", "Stromal II", "Proximal tubule", "Proliferative NPCs", "Endothelial", "Proliferative stromal", "Stromal III", "Stromal IV", "Collecting duct/ distal tubule", "Podocytes", "Unidentified", "Stromal V")
names(new.cluster.ids) <- levels(kidney)
kidney <- RenameIdents(kidney, new.cluster.ids)

tiff("umap_labelled2.tiff", width=1200, height= 1000)
DimPlot(kidney, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

tiff("umap_labelled3.tiff", width=1000, height= 800)
DimPlot(kidney, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

count <- c(1817, 1433, 987, 705, 400, 329, 300, 256, 252, 226, 93, 89, 76)
type <- c("NPC", "Stromal I", "Stromal II", "Proximal tubule", "Proliferative NPCs", "Endothelial", "Proliferative stromal", "Stromal III", "Stromal IV", "Collecting duct/ distal tubule", "Podocytes", "Unidentified", "Stromal V")


# Endothelial sub-cluster
endoc <- subset(kidney, idents = c("Endothelial"))

# Find clusters
endoc <- FindNeighbors(endoc, dims = 1:20)
endoc <- FindClusters(endoc, resolution = 0.5)
head(Idents(endoc), 5)

# UMAP analysis
endoc <- RunUMAP(endoc, dims = 1:20)

tiff("umap_subcluster.tiff", width=1000, height= 800)
DimPlot(endoc, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

tiff("vlnplot_subcluster_NPC.tiff", width=800, height= 800)
VlnPlot(endoc, features = c("SIX1", "SIX2", "CITED1", "PAX2", "ITGA8", "SALL1", "MEIS1", "WT1", "PAX8"))
dev.off()

tiff("featureplot_subcluster_NPC.tiff", width=800, height= 800)
FeaturePlot(endoc, features = c("SIX1", "SIX2", "CITED1", "PAX2", "ITGA8", "SALL1", "MEIS1", "WT1", "PAX8"))
dev.off()

tiff("vlnplot_NPC3.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("SIX1", "SIX2", "CITED1", "PAX2"))
dev.off()

tiff("featureplot_NPC3.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("SIX1", "SIX2", "CITED1", "PAX2"))
dev.off()

tiff("vlnplot_NPC4.tiff", width=800, height= 800)
VlnPlot(kidney, features = c("ITGA8", "MEIS1", "WT1", "PAX8"))
dev.off()

tiff("featureplot_NPC4.tiff", width=800, height= 800)
FeaturePlot(kidney, features = c("ITGA8", "MEIS1", "WT1", "PAX8"))
dev.off()

# Stromal sub-cluster
stromalc <- subset(kidney, idents = c("Stromal I", "Stromal II", "Proliferative stromal", "Stromal III", "Stromal IV", "Stromal V"))

# Find clusters
stromalc <- FindNeighbors(stromalc, dims = 1:20)
stromalc <- FindClusters(stromalc, resolution = 0.5)
head(Idents(stromalc), 5)

# UMAP analysis
stromalc <- RunUMAP(stromalc, dims = 1:20)

tiff("umap_stromalc.tiff", width=1000, height= 800)
DimPlot(stromalc, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

# Nephron and UB only sub-cluster
nephronc <- subset(kidney, idents = c("NPC", "Proximal tubule", "Proliferative NPCs", "Collecting duct/ distal tubule", "Podocytes"))

# Find clusters
nephronc <- FindNeighbors(nephronc, dims = 1:20)
nephronc <- FindClusters(nephronc, resolution = 0.5)
head(Idents(nephronc), 5)

# UMAP analysis
nephronc <- RunUMAP(nephronc, dims = 1:20)

tiff("umap_nephronc.tiff", width=1000, height= 800)
DimPlot(nephronc, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

tiff("umap_endoc.tiff", width=1000, height= 800)
DimPlot(endoc, reduction = "umap", label = TRUE, pt.size = 1)
dev.off()

tiff("dimplotx.tiff", width=800, height= 600)
DimPlot(kidney, reduction = "pca")
dev.off()

# Run TSNE analysis
kidney2 <- RunTSNE(kidney, dims = 1:20)

tiff("tsne.tiff", width=1000, height= 800)
DimPlot(kidney2, reduction = "tsne", label = TRUE, pt.size = 1)
dev.off()


