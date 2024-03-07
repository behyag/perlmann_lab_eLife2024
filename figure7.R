
suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(future)
  require(scales)
  library(RColorBrewer)
  library("readxl")
  library(dplyr)
  library(dendextend)
  library('conover.test')
  library("devtools")
  library(ggpmisc)
  library(VennDiagram)
  library(randomcoloR)
})


### figure 7 + supplements 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

# remove Louvain resolutions from the merged object. 
sobj$seurat_clusters <- NULL

to.remove <- sapply(grep('snn', colnames(sobj@meta.data), value = T),
                    function(x) c(paste(x, collapse = ",")))

for(i in to.remove) {
  sobj[[i]] <- NULL
}

# split the dataset into a list of two seurat objects based on condition 
list <- SplitObject(sobj, split.by = "condition")

# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = SCTransform, variable.features.n=1000)

features <- SelectIntegrationFeatures(object.list = list, nfeatures = 1000)

list <- PrepSCTIntegration(object.list = list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features)

integObj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(integObj) <- "integrated"

integObj <- RunPCA(integObj, npcs = 100, verbose = FALSE)

ElbowPlot(integObj, reduction = "pca", ndims = 100) + ggtitle("mDA_hvg1k_integrated")

# non-linear dim reduction:

DefaultAssay(integObj) <- "integrated"

integObj <- RunUMAP(integObj, reduction = "pca", dims = 1:30, verbose = TRUE)

integObj <- RunTSNE(integObj, reduction = "pca", dims = 1:30, verbose = TRUE)

## find neighbors and find clusters 
integObj <- FindNeighbors(integObj, reduction = "pca", dims = 1:30)

# graph clustering Seurat 
integObj  <- FindClusters(integObj, resolution = seq(1, 4.5, by = 0.5), n.start = 100, n.iter = 100)


## kmeans clustering 
# get the pca embeddings matrix of the desired pca range: 1:30

pcmat <- Embeddings(integObj, reduction = 'pca')[,1:30] 

set.seed(69)
km71 <- kmeans(pcmat, 71, nstart = 100, iter.max = 1000, algorithm="MacQueen")

## add kmeans$clsuter to Seurat Obj metadata 
integObj@meta.data$kmeans71_integ <- km71$cluster 

## get the centroid of the kmeans for hierarchical clustering plus show labels. 
hc71 <- hclust(dist(km71$centers), method = "ward.D2")


# Manual Annotation of integrated territories 
# based on membership (Frequency) of integrated cluster to merged territories. 
# see supplementary file 7 sheet = mergedTERs_integCLUSTERs
# ambiguous clusters were resolved using integrated dendrogram and markers expression. 

anno.df <- as.data.frame(table(integObj$territory, integObj$kmeans71_integ))

# annotation file uploaded, added to integrated object. 
anno.df <- read_excel("/path/supplementary_file7/integrated_territories.xlsx", sheet = 2)
anno.df <- as.data.frame(anno.df)

IntegClusters <- as.character(anno.df$integratedClusters) 

IntegTerritories <- as.character(anno.df$`assigned territories`)

integObj$integTerritory <- plyr::mapvalues(
  x = integObj$kmeans71_integ, 
  from = IntegClusters,
  to = IntegTerritories)

saveRDS(integObj, "/path/to/dir/mDA_integrated.rds")

# fig7 A & B 
Idents(sobj) <- "condition"

tiff(file = "/path/to/dir/fig7A_B.tiff", 
     units="cm", width=40, height=20, res=300)

DimPlot(sobj, split.by = "condition", 
        cols = c("intact"="#08519C", "lesion"="#A63603")) + coord_fixed()

dev.off()


# fig7C merged territories in integrated object
TER_cols <- c("costumized_colors")

tiff(file = "/path/to/dir/fig7c.tiff", 
     units="cm", width=25, height=25, res=300)

DimPlot(integObj, group.by = "territory", cols = TER_cols, 
        order = c("ML", "Hy_DA"), pt.size = 0.7) + 
  theme(legend.text = element_text(size = 12, face = "bold")) + 
  coord_fixed() + ggtitle("merged Territory in integrated obj")

dev.off()

# fig7 D & E Sox6 & Calb1 in merged object split by condition:

# same scale in feature plot 
sapply(grep("Sox6",rownames(sobj@assays$RNA@data),value = TRUE),
       function(x) max(sobj@assays$RNA@data[x,]))

sapply(grep("Calb1",rownames(sobj@assays$RNA@data),value = TRUE),
       function(x) max(sobj@assays$RNA@data[x,]))

sapply(grep("Sox6",rownames(sobj@assays$RNA@data),value = TRUE),
       function(x) min(sobj@assays$RNA@data[x,]))

sapply(grep("Calb1",rownames(sobj@assays$RNA@data),value = TRUE),
       function(x) min(sobj@assays$RNA@data[x,]))

# fixed scale range based on the max & min values of the two genes from above. 

fix.sc <- scale_color_gradientn( colours = c('grey90', 'blue'),  limits = c(0, 4))

tiff(file = "/path/to/dir/fig7D_E.tiff", 
     units="cm", width=25, height=25, res=300)

p1 <- FeaturePlot(sobj, features = c('Sox6', 'Calb1'), 
                  split.by = 'condition', order = T, combine = FALSE) 

p2 <- lapply(p1, function (x) x + fix.sc) 

CombinePlots(p2) 

dev.off()

# Fig7 F highlight ML clusters from merged object in the integrated data
Idents(integObj) <- "territory"

ML_mergedData <- WhichCells(integObj, idents = 'ML')

tiff(file = "/path/to/dir/fig7F.tiff", 
     units="cm", width=25, height=25, res=300)

DimPlot(integObj, cells.highlight = ML_mergedData) + 
  scale_color_manual(labels = c("others", "ML_mergedData"), 
                     values = c("grey", "#08306B")) + 
  theme(legend.text = element_text(size = 10, face = "bold")) + 
  coord_fixed() + ggtitle("merged ML in integrated")

dev.off()


# fig supplement 6-2 A related to integrated object 
# plot UMAP Integrated cluster labels BUT with Territory colors, based on their membership 

clusterCols <- c("costumized_TERcolors")

p <- DimPlot(integObj, group.by = "kmeans71_integ", cols = clusterCols) + coord_fixed() + NoLegend()

tiff(file = "/path/to/dir/figSuppl6-2A.tiff", 
     units="cm", width=25, height=25, res=300)

LabelClusters(plot = p, id = 'kmeans71_integ', color = 'white', size=4, fontface = 'bold', 
              box = T, max.overlapp = Inf, repel=T) 

dev.off()

# fig supplement 6-2 C  VUL module in integrated clusters 
tiff(file = "/path/to/dir/figSuppl6-2C.tiff", 
     units="cm", width=40, height=10, res=300)

VlnPlot(integObj, features = "newVulmodule1", cols = clusterCols, 
        group.by = "kmeans71_integ", pt.size = 0, combine = TRUE) + 
  stat_summary(fun = mean, geom='point', size = 6, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "vulnerable module integ clusters")

dev.off()


# fig supplement 6-2 D  RES module in integrated clusters 

tiff(file = "/path/to/dir/figSuppl6-2D.tiff", 
     units="cm", width=40, height=10, res=300)

VlnPlot(integObj, features = "newResmodule1", cols = clusterCols, 
        group.by = "kmeans71_integ", pt.size = 0, combine = TRUE) + 
  stat_summary(fun = mean, geom='point', size = 6, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "resilient module integ clusters")

dev.off()

# fig supplement 7-1 B 
# plot ML-enriched markers across territories (supplementary file 8)

markers <- c("Atf3", "Creb5", "Xirp2", "Cd9", "Ecel1", "Clic4", "Sprr1a", 
             "Mmp12", "Hrk", "Akain1", "Cd44", "Ighm", "Pim1", "Tnfrsf12a", 
             "Nms", "Igf2bp2", "Arid5a", "Mustn1", "Cd109", "Lifr", "Qrfpr",
             "Pde8a", "Ell2",  "D5Ertd615e", "Bcl2l11", "Dusp10", "Ngf", 
             "Spata13", "Lamb3", "Fbln5")

tiff(file = "/path/to/dir/figSuppl7-1B.tiff", 
     units="cm", width=20, height=20, res=300)

DotPlot(sobj, assay = "RNA", features = markers, group.by = "territory") + coord_flip() + 
  theme(axis.text.x = element_text(size = 14, face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 14, face = "bold"))

dev.off()


# fig supplement 7-1 C
# visualize the overlap of gene sets with VennDiagram package

# import the markers:

# lesion vs intact DE 
lvi <- read.csv("/path/to/dir/lesion_intact_markers.csv", 
                        header = TRUE, sep = ",", quote = "\"")

# ML clusters markers 
ml_up <- read.csv("/path/to/dir/ML_up.csv", 
                  header = TRUE, sep = ",", quote = "\"")

ml_down <- read.csv("/path/to/dir/ML_down.csv", 
                    header = TRUE, sep = ",", quote = "\"")

# get the gene vectors:
lvi_up <- lvi[lvi$p_val_adj < 0.05 & lvi$avg_log2FC > 0, ]$X
lvi_down <- lvi[lvi$p_val_adj < 0.05 & lvi$avg_log2FC < 0, ]$X
ml_up <- ml_up[ml_up$p_val_adj < 0.05, ]$X
ml_down <- ml_down[ml_down$p_val_adj < 0.05, ]$X

# set the colors:
brewer.pal(3, 'Paired')
# "#A6CEE3" "#1F78B4" "#B2DF8A"
show_col(brewer.pal(3, 'Paired'))
myCol <- c("#A6CEE3", "#B2DF8A")

venn.diagram(
  x = list(ml_up, lvi_up),
  category.names = c("ML upregulated" , "lesioned upregulated"),
  filename = '/path/to/dir/figSupp7-1C-up.tiff',
  output=TRUE,
  imagetype="tiff" ,
  height = 700 , 
  width = 700 , 
  resolution = 2000,
  compression = "lzw",
  main = 'up-regulayed genes',
  main.fontface = 'bold', 
  main.cex = 0.05,
  lwd = 0,
  lty = 'blank', 
  fill = myCol,
  cex = 0.25,
  fontface = "bold",
  ext.text = TRUE,
  ext.line.lwd = 0.1, 
  ext.dist = -0.09,
  ext.length = 0.8,
  cat.cex = 0.05,
  cat.fontface = "bold", 
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.2), 
  margin = 0.05)


venn.diagram(
  x = list(ml_down, lvi_down),
  category.names = c("ML down-regulated" , "lesion down-regulated"),
  filename = '/path/to/dir/figSupp7-1C-down.tiff',
  output=TRUE,
  imagetype="tiff" ,
  height = 700 , 
  width = 700 , 
  resolution = 2000,
  compression = "lzw",
  main = 'down-regulayed genes',
  main.fontface = 'bold', 
  main.cex = 0.05,
  lwd = 0,
  lty = 'blank', 
  fill = myCol,
  cex = 0.25,
  fontface = "bold",
  ext.text = 1,
  ext.line.lwd = 0.1, 
  ext.dist = 0.1,
  ext.length = 0.8,
  cat.cex = 0.05,
  cat.fontface = "bold", 
  cat.default.pos = "outer",
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.2), 
  margin = 0.05)

# fig supplement 7-1 D & E 
# mostly_lesioned (ML) clusters enriched markers GSE analysis
# only genes unique to ML clusters and NOT the ones also enriched in lesion_intact DE: 

ml_up_unique <- read.csv("/path/to/dir/ML_up_unique.csv", 
                         header = TRUE, sep = ",", quote = "\"")

ml_down_unique <- read.csv("/path/to/dir/ml_down_unique.csv", 
                           header = TRUE, sep = ",", quote = "\"")

# get the gene names column:  
ml_up_unique <- ml_up_unique$x
ml_down_unique <- ml_down_unique$x

# next,  Enrichr for analysis 
install.packages("enrichR")

library(enrichR)

# list Enrichr sites:

listEnrichrSites()

dbs <- listEnrichrDbs()

websiteLive <- TRUE

if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# select the desired gene_set_libraries from the Enrichr databases:  

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")

if (websiteLive) {
  enriched_up <- enrichr(ml_up_unique, dbs)
}

if (websiteLive) {
  enriched_down <- enrichr(ml_down_unique, dbs)
}

# you may make and write the results tables.

Up_mol <- if (websiteLive) enriched_up[["GO_Molecular_Function_2021"]]

Up_cell <- if (websiteLive) enriched_up[["GO_Cellular_Component_2021"]]

Up_bio <- if (websiteLive) enriched_up[["GO_Biological_Process_2021"]]

write.csv(Up_mol, '/path/to/dir/up_GO_Molecular_Function_2021.csv')

write.csv(Up_cell, '/path/to/dir/up_GO_Cellular_Component_2021.csv')

write.csv(Up_bio, '/path/to/dir/up_GO_Biological_Process_2021.csv')

down_mol <- if (websiteLive) enriched_down[["GO_Molecular_Function_2021"]]

down_cell <- if (websiteLive) enriched_down[["GO_Cellular_Component_2021"]]

down_bio <- if (websiteLive) enriched_down[["GO_Biological_Process_2021"]]

write.csv(down_mol, '/path/to/dir/down_GO_Molecular_Function_2021.csv')
write.csv(down_cell, '/path/to/dir/down_GO_Cellular_Component_2021.csv')
write.csv(down_bio, '/path/to/dir/down_GO_Biological_Process_2021.csv')

## Plot Enrichr GO-BP output. (Plotting function contributed by I-Hsuan Lin)

p1 <- if (websiteLive) plotEnrich(enriched_up[['GO_Biological_Process_2021']], 
                                  showTerms = 20, numChar = 90, y = "Count", orderBy = "P.value",
                                  title = 'up-regulated GO_Biological_Process_2021')

p1 + theme(plot.title = element_text(size = 13, face = "bold"), 
           axis.text.x = element_text(size = 10, face = "bold"), 
           axis.text.y = element_text(size = 10, face = "bold"))


p2 <- if (websiteLive) plotEnrich(enriched_down[['GO_Biological_Process_2021']], 
                                  showTerms = 20, numChar = 90, y = "Count", orderBy = "P.value",
                                  title = 'down-regulated GO_Biological_Process_2021')

p2 + theme(plot.title = element_text(size = 13, face = "bold"), 
           axis.text.x = element_text(size = 10, face = "bold"), 
           axis.text.y = element_text(size = 10, face = "bold"))



