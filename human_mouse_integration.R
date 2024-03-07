
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(stringr)
  library(sctransform)
  library(future)
  require(scales)
  library(RColorBrewer)
  library("readxl")
  library(dplyr)
  library(dendextend)
  library(patchwork)
  library(ggplot2)
  library(Orthology.eg.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})

###  re-analysis of human data Kamath et al. NatureNeuroscience 2022

### Generation of human dataset

# data (count matrix, barcodes and features) downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265
# GSE178265_Homo_bcd.tsv.gz	
# GSE178265_Homo_features.tsv.gz	
# GSE178265_Homo_matrix.mtx.gz

# data files renamed and seurat object made. 

km <- Read10X("/path/to/files/data")

km <- CreateSeuratObject(km, project = 'km')

# metadata_PD tsv file downloaded from:
# https://singlecell.broadinstitute.org/single_cell/study/SCP1768/single-cell-genomic-profiling-of-human-dopamine-neurons-identifies-a-population-that-selectively-degenerates-in-parkinsons-disease-single-nuclei-data#study-download

# this is the metadata of all nuclei from all species in this study. 

metadata <- read.table(file = '/path/to/dir/METADATA_PD.tsv', sep = '\t', header = T)

# get cell barcodes from the km Seurat object and subset the metadata based on that:

barcodes <- FetchData(km, vars = 'ident')

# turn barcodes' row names into a vector:
cellbc <- row.names(barcodes)

# subset metadata based on this vector:

metadf <- metadata[metadata$NAME %in% cellbc, ] 

write.table(metadf, file="/path/to/dir/Metadata_human.tsv", 
            quote=FALSE, sep='\t', col.names = TRUE)

# it's only human now:
unique(metadf$species__ontology_label )
#  [1] "Homo sapiens"

# set the barcodes (NAME column) as the row names of the metadf: 
row.names(metadf) <- metadf$NAME 

km <- AddMetaData(object = km, metadata = metadf, col.name = NULL)

km <- PercentageFeatureSet(km, "^MT-", col.name = "percent_mito", assay = "RNA")

km <- PercentageFeatureSet(km, "^RP[SL]", col.name = "percent_ribo", assay = "RNA")

km <- CellCycleScoring(km, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

km$UMIperGene <- km$nCount_RNA / km$nFeature_RNA

### gene & cell filtering 

selected_cells <- WhichCells(km, expression = nCount_RNA >= 650 & percent_mito <= 10 & UMIperGene >= 1.2 )

selected_genes <- rownames(km)[Matrix::rowSums(km) > 1]

km <- subset(km, features = selected_genes, cells = selected_cells)

# remove MALAT1

km <- km[ ! grepl("MALAT1", rownames(km)), ]

saveRDS(km, "/path/to/dir/km.rds")

### subset dataset based on condition & organ 

table(km$disease__ontology_label )   
#  Lewy body dementia       normal      Parkinson disease 
#              67466        231530                135344 

# first, subset normal:

km <- SetIdent(km, value = 'disease__ontology_label')

kmsub1 <- subset(km, idents = 'normal')

kmsub1 <- SetIdent(kmsub1, value = 'organ__ontology_label')

kmsub1 <- subset(kmsub1, idents = 'substantia nigra pars compacta')

saveRDS(kmsub1, "/path/to/dir/km_cntl.rds")

### subset diseased individuals 

kmsub2 <- subset(km, idents = c('Parkinson disease', 'Lewy body dementia'))

kmsub2 <- SetIdent(kmsub2, value = 'organ__ontology_label')

kmsub2 <- subset(kmsub2, idents = 'substantia nigra pars compacta')

saveRDS(kmsub2, "/path/to/dir/km_PDLBD.rds")


### integration of control SNpc dataset 

# load all data from previous step 
mydata <- readRDS("/path/to/dir/km_cntl.rds")

# split the dataset by donor_id
my.list <- SplitObject(mydata, split.by = 'donor_id')
gc()
# normalize and identify variable features for each dataset independently
my.list <- lapply(X = my.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

gc()

# select features that are repeatedly variable across datasets for integration nfeatures = 2000 default 
features <- SelectIntegrationFeatures(object.list = my.list, nfeatures = 2000)

gc()

# scale and run PCA on each object in the list. 

my.list <- lapply(X = my.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

gc()

# Perform integration
my.anchors <- FindIntegrationAnchors(object.list = my.list, reduction = "rpca", dims = 1:50)

rm(my.list)

gc()

# creates an 'integrated' data assay 
mydata <- IntegrateData(anchorset = my.anchors, dims = 1:50)

gc()

# Perform integrated analysis
DefaultAssay(mydata) <- "integrated"

# Run the standard workflow for visualization and clustering
mydata <- ScaleData(mydata, verbose = FALSE)

mydata <- RunPCA(mydata, npcs = 50, verbose = FALSE)

tiff(file = "/path/to/dir/elbowplot.tiff", 
     units="cm", width = 50, height = 50, res = 300)

ElbowPlot(mydata, reduction = "pca", ndims = 50) + ggtitle("km_Cntl integrated")

dev.off()

gc()

mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:30)

mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30)

mydata <- FindClusters(mydata, resolution = seq(0.1, 1.5, by=0.1), n.start = 100, n.iter = 100)

saveRDS(mydata, "/path/to/dir/km_cntl_integ_snpc.rds")

print(sessionInfo())


### integration of PD LBD dataset 

# load all data from previous step 
mydata <- readRDS("/path/to/dir/km_PDLBD.rds")

print(mydata)

# split the dataset by donor_id
my.list <- SplitObject(mydata, split.by = 'donor_id')
gc()
# normalize and identify variable features for each dataset independently
my.list <- lapply(X = my.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

gc()

# select features that are repeatedly variable across datasets for integration nfeatures = 2000 default 
features <- SelectIntegrationFeatures(object.list = my.list, nfeatures = 2000)

gc()

# scale and run PCA on each object in the list. 

my.list <- lapply(X = my.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

gc()

### Perform integration
my.anchors <- FindIntegrationAnchors(object.list = my.list, reduction = "rpca", dims = 1:50)

rm(my.list)

gc()

#  creates an 'integrated' data assay 
mydata <- IntegrateData(anchorset = my.anchors, dims = 1:50)

gc()

# Perform integrated analysis
DefaultAssay(mydata) <- "integrated"

# Run the standard workflow for visualization and clustering
mydata <- ScaleData(mydata, verbose = FALSE)

mydata <- RunPCA(mydata, npcs = 50, verbose = FALSE)

tiff(file = "/path/to/dir/elbowplot.tiff", 
     units="cm", width = 50, height = 50, res = 300)

ElbowPlot(mydata, reduction = "pca", ndims = 50) + ggtitle("km_PDLBD integrated")

dev.off()

gc()

mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:30)

mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30)

mydata <- FindClusters(mydata, resolution = seq(0.1, 1.5, by=0.1), n.start = 100, n.iter = 100)

saveRDS(mydata, "/path/to/dir/km_PDLBD_integrated.rds")

print(sessionInfo())


###  integration of control with PDLBD 
# default assay was set to integrated for both control and PDLBD objects, 
# to enable use of the the already calculated PCA in each. 
# Reciprocal PCA (RPCA) was used. 

cnt <- readRDS("/path/to/dir/km_cntl_integ_snpc.rds")

DefaultAssay(cnt) <- 'integrated'

pd <- readRDS("/path/to/dir/km_PDLBD_integrated.rds")

DefaultAssay(pd) <- 'integrated'

# Perform integration

my.anchors <- FindIntegrationAnchors(object.list = list(cnt, pd), reduction = "rpca", dims = 1:50)

rm(cnt)
rm(pd)
gc()

#  create an 'integrated' data assay
mydata <- IntegrateData(anchorset = my.anchors)

gc()

# Perform integrated analysis
DefaultAssay(mydata) <- "integrated"

# Run the standard workflow for visualization and clustering
mydata <- ScaleData(mydata, verbose = FALSE)

mydata <- RunPCA(mydata, npcs = 50, verbose = FALSE)

tiff(file = "/path/to/dir/elbowplot.tiff", 
     units="cm", width = 40, height = 40, res = 300)

ElbowPlot(mydata, reduction = "pca", ndims = 50) + ggtitle("All integrated")

dev.off()

gc()

mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:30)

mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30)

mydata <- FindClusters(mydata, resolution = seq(0.1, 1.5, by=0.1), n.start = 100, n.iter = 100)

saveRDS(mydata, "/path/to/dir/km_All_integrated.rds")

### subset DA clusters based on canonical dopaminergic markers expression 
# clusters #4 and #10 at integrated_snn_res.0.1 resolution are dopaminergic

mydata <- SetIdent(mydata, value = 'integrated_snn_res.0.1')

km_DA <- subset(mydata, idents = c(4, 10))

table(km_DA$Status )

#  Ctrl   LBD    PD 
# 17039  4642  3322 

saveRDS(km_DA, "/path/to/dir/kmAll_DA.rds")

print(sessionInfo())


####### Integration of human and mouse dataset

# counts matrix of 'RNA' assay from human DA dataset
km_DA <- readRDS("/path/to/dir/kmAll_DA.rds")

DefaultAssay(km_DA) <- 'RNA'

Hscounts <- GetAssayData(km_DA, slot = 'counts')

write.table(Hscounts, file= '/path/to/dir/Hs_counts.tsv', sep='\t', col.names = T)

# gene map function 
mapfun <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  out
}

# keys = genes symbols 
z <- keys(org.Mm.eg.db, "SYMBOL")

Hs_Mm_mappedgenes <- mapfun(z)

write.csv(Hs_Mm_mappedgenes, '/path/to/dir/Hs_Mm_mappedgenes.csv', row.names = TRUE)

# create the human count matrix, with Orthologous mouse gene IDs
mmhs <- read.csv('/path/to/dir/Hs_Mm_mappedgenes.csv', sep = ',')

#  Row names in the metadata need to match the column names of the counts matrix, header = T

Hscounts <- read.table(file = '/path/to/dir/Hs_counts.tsv', sep = '\t', header = T) 

jgenes <- intersect(mmhs$Human_symbol, rownames(Hscounts))

# two subsets, genes with / without matching mouse symbol 
mdf <- Hscounts[rownames(Hscounts) %in% jgenes, ]

mdf2 <- Hscounts[!(rownames(Hscounts) %in% jgenes), ]

mmgenes <- mmhs[mmhs$Human_symbol %in% jgenes, ]

# re-arrange columns 
mmgenes <- mmgenes[, c(2, 5)]

# re-order mmgenes based on mdf rows:
mmgenes <- mmgenes[match(rownames(mdf), mmgenes$Human_symbol ),]

identical(rownames(mdf), mmgenes$Human_symbol )
# [1] TRUE

mdf <- cbind(mdf, mmgenes)

rownames(mdf) <- mdf$Mouse_symbol

mdf <- mdf[, -c(25004:25005)] 

# now the second part of the count matrix df: mdf2

# create a new column based on row names to modify them and reuse as row names 
mdf2$c1 <- rownames(mdf2) 

identical(rownames(mdf2), mdf2$c1 )
#[1] TRUE

mdf2$c1 <- tolower(mdf2$c1 ) 

# First letter to upper case:
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mdf2$c1 <- firstup(mdf2$c1 ) 

mdf2$c1 <- make.unique(mdf2$c1 ) 

mdf2$c1  <-  gsub('\\.', '-', mdf2$c1 )

rownames(mdf2) <- NULL

rownames(mdf2) <- mdf2$c1 

mdf2$c1 <- NULL 

# join the two data frames to create the whole count matrix based on mouse orthologs and unique human genes

merged_mdf <- rbind.data.frame(mdf, mdf2)

dim(merged_mdf)
#[1] 38945 25003

#same dim as the original count matrix
dim(Hscounts)
# [1] 38945 25003

write.table(merged_mdf, file= '/path/to/dir/Hs_counts_Mm.tsv', sep='\t', col.names = T)

# pull the metadata info from Seurat object 
metadata <- km_DA@meta.data 
dim(metadata)
#  [1] 25003    44

colnames(metadata)

# remove extra metdata columns 
# remove all "integrated_snn_res...." columns:
to.go <- sapply(grep('snn', colnames(metadata), value = T),
                function(x) c(paste(x, collapse = ",")))

for(i in to.go) {
  metadata[[i]] <- NULL
}

# remove "Name" and other unwanted columns:
to.go <- c('seurat_clusters', 'NAME', 'libname', 'species', 'organ', 'library_preparation_protocol', 
           'orig.ident', 'disease' )

for(i in to.go) {
  metadata[[i]] <- NULL
}

# set donor_ID as orig.ident:
metadata$orig.ident <- metadata$donor_id  

identical(metadata$donor_id, metadata$orig.ident  )
# [1] TRUE

metadata$donor_id <- NULL 

metadata <- metadata[, c(21, 1:20)] 

colnames(metadata)

# re-arrange metadata to make it more compatible with mouse metadata and vice versa 

metadata$species <- 'Homo_sapiens' 

metadata$species__ontology_label <- NULL

metadata <- metadata[, c(1:3, 21, 4:20)] 

metadata$disease__ontology_label <- NULL

metadata$lib_prep_protocol <- metadata$library_preparation_protocol__ontology_label 

metadata$library_preparation_protocol__ontology_label <- NULL 

metadata <- metadata[, c(1:4, 7, 11, 14:20, 5, 6, 8:10, 12, 13)] 

colnames(metadata)

# load the counts.tsv file, generated above:  

counts <- read.table(file = '/path/to/dir/Hs_counts_Mm.tsv', sep = '\t', header = T) 

# in the counts df, from which matrix.mtx is generated, cell IDs have a <.> but, 
# in metadata df, cell IDs have <-> instead. so for compatibility, <.>  must be changed to <->

colnames(counts) <-  gsub('\\.', '-', colnames(counts) )

identical(rownames(metadata), colnames(counts))
#  [1] TRUE

# re-write the new, updated counts table 
write.table(counts, file= '/path/to/dir/HsMmcounts.tsv', sep='\t', col.names = T)

# generate the 3 files required for Read10x() (in a new directory) 
write(x = rownames(counts), file = '/path/to/new_dir/features.tsv', sep = '\t')

write(x = colnames(counts), file = '/path/to/new_dir/barcodes.tsv', sep = '\t')

mat <- data.matrix(counts)
sp.mat <- Matrix(mat, sparse = T)
writeMM(obj = sp.mat, file = '/path/to/new_dir/matrix.mtx')

# in the written features.tsv file, there's only one column with gene names
#  but read10x(), takes column #2 for gene names by default, so gene.column = 1          

HsDA <- Read10X("/path/to/new_dir", gene.column = 1, cell.column = 1)

# Create Seurat Object for human dopaminergic data
HsDA <- CreateSeuratObject(HsDA, meta.data = metadata, project = 'HsDA')

View(HsDA@meta.data)

saveRDS(HsDA, '/path/to/dir/HsDA.rds')

### load the mouse mDA dataset generated in mDA.R 

mda <- readRDS("/path/to/dir/mDA.rds")

DefaultAssay(mda) <- 'RNA'

# get the counts matrix of 'RNA' assay from mouse mDA dataset

Mmcounts <- GetAssayData(mda, slot = 'counts')

write.table(Mmcounts, file= '/path/to/dir/Mm_counts.tsv', sep='\t', col.names = F)

#  change sparse matrix to dataframe

Mmcounts <- as.data.frame(as.matrix(Mmcounts))

dim(Mmcounts)
# [1] 26497 33052

# pull the metadata info from Seurat object 
metadata <- mda@meta.data 
dim(metadata)
#  [1] 33052    64

colnames(metadata)

# function to remove all "integrated_snn_res...." columns:

to.go <- sapply(grep('snn', colnames(metadata), value = T),
                function(x) c(paste(x, collapse = ",")))

for(i in to.go) {
  metadata[[i]] <- NULL
}

to.go <- sapply(grep('kmeans', colnames(metadata), value = T),
                function(x) c(paste(x, collapse = ",")))

for(i in to.go) {
  metadata[[i]] <- NULL
}

metadata <- metadata[, -c(13:30)] 

colnames(metadata)

dim(metadata)
# 33052    16

### Create Seurat Object from the counts df and the metadata df: 

# Row names in the metadata need to match the column names of the counts matrix

identical(rownames(metadata), colnames(Mmcounts))
#  [1] TRUE

MmDA <- CreateSeuratObject(counts = Mmcounts, meta.data = metadata, project = 'MmDA')

View(MmDA@meta.data )

# re-arrange metadata to make it more compatible for integration 

MmDA$species <- 'Mus_musculus' 

MmDA@meta.data <- MmDA@meta.data[, c(1:3, 17, 4:16)]

MmDA$UMIperGene <- MmDA$UMIsPerGene 

MmDA$UMIsPerGene <- NULL

MmDA$Status <- MmDA$condition  

MmDA$condition <- NULL

MmDA$GenesPerUMI <- NULL 

MmDA$sex <- 'female'

MmDA@meta.data <- MmDA@meta.data[, c(1:4, 17, 16, 6:10, 15, 5, 11:14)]  

names(MmDA@meta.data )

saveRDS(MmDA, '/path/to/dir/MmDA.rds')


### integration 

# load human & mouse datasets generated above
MmDA <- readRDS('/path/to/dir/MmDA.rds')

HsDA <- readRDS('/path/to/dir/HsDA.rds')

mylist <- list(HsDA, MmDA)

# normalize and identify variable features for each dataset independently
mylist <- lapply(X = mylist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = mylist)

myanchors <- FindIntegrationAnchors(object.list = mylist, anchor.features = features)

rm(MmDA)
rm(HsDA)
gc()

#  create the 'integrated' data assay
mydata <- IntegrateData(anchorset = myanchors)

# Perform integrated analysis
DefaultAssay(mydata) <- "integrated"

# Run the standard workflow for visualization and clustering
mydata <- ScaleData(mydata, verbose = FALSE)

mydata <- RunPCA(mydata, npcs = 50, verbose = FALSE)

tiff(file = "/path/to/dir/elbowplot.tiff", 
     units="cm", width = 30, height = 30, res = 300)

ElbowPlot(mydata, reduction = "pca", ndims = 50) + ggtitle("HsMm integrated")

dev.off()

mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:30)

mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:30)

mydata <- FindClusters(mydata, resolution = seq(0.1, 1.5, by=0.1), n.start = 100, n.iter = 100)

## create a new metadata column for both status and sample ID
mydata <- SetIdent(mydata, value = 'Status')

mydata$status.sample <- paste(Idents(mydata), mydata$orig.ident, sep = "_")

saveRDS(mydata, "/path/to/dir/DAHsMmIntegrated.rds")

### figure 8 

ord <- c("Ctrl", "PD", "LBD", "intact", "lesion")

mydata$Status <- factor(mydata$Status, levels = ord ) 

p <-  DimPlot(mydata, split.by = 'Status', group.by = 'Status', ncol = 3) + coord_fixed()

tiff(file = "/path/to/dir/fig8.tiff", 
     units="cm", width=50, height=25, res=300)

p + theme(legend.text = element_text(size = 18, face = 'bold'), 
          strip.text.x = element_text(size = 18, face = "bold"))

dev.off()

DimPlot(mydata, group.by = 'Status') + coord_fixed()

print(sessionInfo())









