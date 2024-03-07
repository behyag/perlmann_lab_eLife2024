
suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(RColorBrewer)
})


### figure Supplement 8-1 

# Mapping and annotating human query dataset with mouse as ref dataset

# ref: mouse intact-only mDA territories
# query: human control mDA from SNpc

# load the two datasets 
MmDA <- readRDS('/path/to/dir/MmDA.rds')
HsDA <- readRDS('/path/to/dir/HsDA.rds')

HsDA <-SetIdent(HsDA, value = 'Status')
MmDA <-SetIdent(MmDA, value = 'Status')

# subset control from human and intact from mouse datasets. 
HsDA <- subset(HsDA, idents = 'Ctrl')
MmDA <- subset(MmDA, idents = 'intact')

HsDA <- NormalizeData(HsDA, normalization.method = "LogNormalize", scale.factor = 10000)

HsDA <- ScaleData(HsDA)

# for mouse data, rename territories & unassigned clusters 
# based on being dopaminergic or non-dopaminergic (nonmDA) 

# subset only the dopaminergic ones 

MmDA <-SetIdent(MmDA, value = 'territory')

unique(levels(MmDA$territory))

MmDA$TER <- plyr::mapvalues(
  x = MmDA$territory, 
  from = c("ML", "Gad2", "Fbn2", "Pcsk6", "Pdia5", "Ebf1", "Otx2", "mODC", "Sox6", "Lef1", 
           "HPT", "Vglut1", "36", "35", "STN_Vglut2", "GABA", "51" ),
  to = c("ML", "Gad2", "Fbn2", "Pcsk6", "Pdia5", "Ebf1", "Otx2", "nonmDA", "Sox6", "nonmDA", 
         "nonmDA", "nonmDA", "nonmDA", "nonmDA", "nonmDA", "nonmDA", "nonmDA"))

MmDA <-SetIdent(MmDA, value = 'TER')

MmDA <- subset(MmDA, idents = 'nonmDA', invert=T)

MmDA <- NormalizeData(MmDA, normalization.method = "LogNormalize", scale.factor = 10000)

MmDA <- FindVariableFeatures(MmDA, selection.method = "vst", nfeatures = 2000)

MmDA <- ScaleData(MmDA)

MmDA <- RunPCA(MmDA, npcs = 50, verbose = F)

ElbowPlot(MmDA, reduction = "pca", ndims = 50) + ggtitle("MmDA intact non-mDA TERs removed")

### Unimodal UMAP Projection
#   enable projection of a query onto the reference UMAP structure. 

MmDA <- RunUMAP(MmDA, dims = 1:30, reduction = "pca", return.model = TRUE)

myanchors <- FindTransferAnchors(reference = MmDA, query = HsDA, 
                                 dims = 1:30, k.anchor = 5, 
                                 k.filter = NA, 
                                 reference.reduction = "pca")

HsDA <- MapQuery(anchorset = myanchors, reference = MmDA, query = HsDA,
                 refdata = MmDA$TER , reference.reduction = "pca", 
                 reduction.model = "umap")

tiff(file = "/path/to/dir/figS8a.tiff", 
     units="cm", width=50, height=25, res=300)

p1 <- DimPlot(MmDA, reduction = "umap", group.by = "TER", label = F) + 
  coord_fixed() + ggtitle("Reference territory annotations") + 
  theme(legend.text = element_text(size = 18, face = 'bold')) + 
  guides(colour = guide_legend(override.aes = list(size=9, alpha = 1)))

p1 <- LabelClusters(p1, id = 'TER', fontface = 'bold', size = 6, repel = T)
p1
dev.off()

tiff(file = "/path/to/dir/figS8b.tiff", 
     units="cm", width=50, height=25, res=300)

p2 <- DimPlot(HsDA, reduction = "ref.umap", group.by = "predicted.id", label = F, order = 'Sox6') + 
  coord_fixed() + ggtitle("Query transferred labels") + 
  theme(legend.text = element_text(size = 18, face = 'bold')) + 
  guides(colour = guide_legend(override.aes = list(size=9, alpha = 1)))

p2 <- LabelClusters(p2, id = 'predicted.id', fontface = 'bold', size = 6, repel = T)

p2
dev.off()

HsDA <- SetIdent(HsDA, value = 'predicted.id')

tiff(file = "/path/to/dir/figS8.tiff", 
     units="cm", width=25, height=25, res=300)

DimPlot(HsDA, reduction = "ref.umap", split.by = "predicted.id", label = F,
        cols = c('Sox6'='#FF61CC', 'Gad2'='#CD9600', 'Otx2'='#C77CFF', 'Pdia5'='#00BFC4'), ncol = 2) + 
  coord_fixed() + ggtitle("Query transferred labels") + 
  theme(legend.text = element_text(size = 18, face = 'bold'), 
        strip.text.x = element_text(size = 18, face = "bold") ) + 
  guides(colour = guide_legend(override.aes = list(size=9, alpha = 1)))

dev.off()

# figure supplement 8-1 D pie chart of human query nuclei with mouse label transferred 

tb1 <- table(HsDA$predicted.id )

lbls <- paste(names(tb1), "\n", tb1, sep="")

tiff(file = "/path/to/dir/figS8D.tiff", 
     units="cm", width=25, height=25, res=300)

pie(tb1, labels = lbls, 
    col = c('Gad2'='#CD9600', 'Otx2'='#C77CFF', 'Pdia5'='#00BFC4', 'Sox6'='#FF61CC'),
    main="Pie Chart of Human control\n (per mouse mDA territory sizes)", cex = 2) 

dev.off()


## subset only sox6 labeled human nuclei as query, and mouse Sox6 neighborhoods as reference. 

# ref: mouse intact only Sox6 neighborhoods

# query: human control only Sox6-TER labeled nuclei 

# 1. human query: Sox6 labeled transferred
HsDA <-SetIdent(HsDA, value = 'predicted.id')
hssox6 <- subset(HsDA, idents = 'Sox6')
hssox6[['prediction.score.id']] <- NULL
hssox6[['ref.pca']] <- NULL
hssox6[['ref.umap']] <- NULL

# 2. mouse ref: Sox6 territory
MmDA <-SetIdent(MmDA, value = 'TER')
mmsox6 <- subset(MmDA, idents = 'Sox6')

mmsox6 <- NormalizeData(mmsox6, normalization.method = "LogNormalize", scale.factor = 10000)

mmsox6 <- FindVariableFeatures(mmsox6, selection.method = "vst", nfeatures = 2000)

mmsox6 <- ScaleData(mmsox6)

mmsox6 <- RunPCA(mmsox6, npcs = 50, verbose = F)

###   Unimodal UMAP Projection

mmsox6 <- RunUMAP(mmsox6, dims = 1:30, reduction = "pca", return.model = TRUE)

myanchors <- FindTransferAnchors(reference = mmsox6, query = hssox6, dims = 1:30, 
                                 reference.reduction = "pca", k.filter = NA)

hssox6 <- MapQuery(anchorset = myanchors, reference = mmsox6, query = hssox6,
                   refdata = mmsox6$neighborhood, 
                   reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(mmsox6, reduction = "umap", group.by = "neighborhood", label = TRUE, label.size = 5,
              repel = TRUE) + coord_fixed() + 
  ggtitle("Reference neighborhood annotations")

p2 <- DimPlot(hssox6, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 5, repel = T) + coord_fixed() + 
  ggtitle("Query transferred labels")

p1 + p2


hssox6 <- SetIdent(hssox6, value = 'predicted.id')

DimPlot(hssox6, reduction = "ref.umap", split.by = "predicted.id", label = F,
        cols = c('Sox6_NH1'='#F8766D', 'Sox6_NH2'='#C77CFF', 'Sox6_NH4'='#7CAE00', 'Sox6_NH3'='#00BFC4'), 
        ncol = 2) + coord_fixed() + ggtitle("Query transferred labels") +
  theme(legend.text = element_text(size = 18, face = 'bold'), 
        strip.text.x = element_text(size = 18, face = "bold") ) + 
  guides(colour = guide_legend(override.aes = list(size=9, alpha = 1)))


mmsox6 <- SetIdent(mmsox6, value = 'neighborhood')

DimPlot(mmsox6, reduction = "umap", split.by = "neighborhood", label = F, label.size = 5,
        cols = c('Sox6_NH1'='#F8766D', 'Sox6_NH2'='#C77CFF', 'Sox6_NH4'='#7CAE00', 'Sox6_NH3'='#00BFC4'), 
        ncol = 2) + coord_fixed() + ggtitle("Ref Sox6_NH")

# Pie Chart 
tb1 <- table(hssox6$predicted.id )

lbls <- paste(names(tb1), "\n", tb1, sep="")

pie(tb1, labels = lbls,
    main="Pie Chart of Human control\n (per mouse mDA territory sizes)")


print(sessionInfo())

