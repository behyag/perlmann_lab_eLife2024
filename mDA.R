
suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(sctransform)
  library(ggplot2)
  library(ggpubr)
  library(future)
  require(scales)
  library(RColorBrewer)
  library("readxl")
  library(dplyr)
  library(dendextend)
  library('conover.test')
})


# import the all_nuclei dataset, created in allnuclei.R
#  sobj = Seurat Object 
sobj <- readRDS("/path/to/dir/allnuc.rds")

# filter all_nuclei datset based on Th / Slc6a3 (DAT) expression. 

sobj <- subset(sobj, subset = Slc6a3 > 0 | Th > 0)

# remove SCT assay and resolutions left from the previous object (all_nuclei)

sobj[["SCT"]] <- NULL
sobj$seurat_clusters <- NULL

to.remove <- sapply(grep('snn', colnames(sobj@meta.data), value = T),
                function(x) c(paste(x, collapse = ",")))

for(i in to.remove) {
  sobj[[i]] <- NULL
}

# a modeling framework for the normalization and variance stabilization of molecular count data.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform.
sobj <- SCTransform(sobj, variable.features.n = 1000)

seed.use = # set to the random default seed for each Seurat function, unless stated otherwise! 
  
sobj <- RunPCA(sobj, npcs = 50, verbose = FALSE)
ElbowPlot(asobjn, reduction = "pca", ndims = 50)

sobj <- RunTSNE(sobj, dims = 1:30)
sobj <- RunUMAP(sobj, dims = 1:30)

# SNN Graph Construction
sobj <- FindNeighbors(sobj, reduction = "pca", k.param = 20, dims = 1:30)

# k-means clustering 
pcmat <- Embeddings(sobj, reduction = 'pca')[,1:30]

set.seed(69)
km71 <- kmeans(pcmat, 71, nstart = 50, iter.max = 1000, algorithm="MacQueen")

hc71 <- hclust(dist(km71$centers ), method = "ward.D2")

plot(hc71, hang = -1)

sobj@meta.data$kmeans71 <- km71$cluster 

# log-normalization of genes in "RNA" assay 

DefaultAssay(sobj) <- "RNA"

sobj <- NormalizeData(sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

### create a new metadata column for territory x condition 
Idents(sobj) <- 'territory'

sobj$territory.condition <- paste(Idents(sobj), sobj$condition, sep = '_') 

### create a new metadata column for neighborhood x condition 
Idents(sobj) <- 'neighborhood'

sobj$neighborhood.condition <- paste(Idents(sobj), sobj$condition, sep = '_') 

# save the object 
saveRDS(sobj, "/path/to/dir/mDA.rds") 

# cluster markers identification 
sobj <- SetIdent(sobj, value = "kmeans71")

print(length(levels(sobj@active.ident)))

plan("multicore", workers = 12)

options(future.globals.maxSize = 96000 * 1024^2)

for (i in 1:71){
  DE <- assign(paste0(i, "kmeans71"), 
               FindMarkers(sobj, assay = "RNA", slot = "data", ident.1 = i, verbose = FALSE) )
  write.csv(DE, file=paste0("/path/to/dir/c_",i,".csv"))
}


# Territory markers 

DefaultAssay(sobj) <- "RNA"

sobj <- SetIdent(sobj, value = "territory")

plan("multicore", workers = 12)

options(future.globals.maxSize = 96000 * 1024^2)

Markers <- FindAllMarkers(sobj, slot = 'data')

write.csv(Markers, "/path/to/dir/TerritoryMarkers.csv", row.names = TRUE)


# Neighborhood markers 

DefaultAssay(sobj) <- "RNA"

sobj <- SetIdent(sobj, value = "neighborhood")

plan("multicore", workers = 12)

options(future.globals.maxSize = 96000 * 1024^2)

Markers <- FindAllMarkers(sobj, slot = 'data')

write.csv(Markers, "/path/to/dir/NeighborhoodMarkers.csv", row.names = TRUE)


# DE genes between lesion and intact irrespective of clusters. 

DefaultAssay(sobj) <- "RNA"

sobj <- SetIdent(sobj, value = "condition")

Markers <- FindMarkers(sobj, ident.1 = "lesion", ident.2 = "intact", only.pos = FALSE)

write.csv(Markers, "/path/to/dir/lesion_intact_markers.csv", row.names = TRUE)


### cell loss
# identification of sub-clusters for k-means71, as the basis for the calculation of normalized cell loss. 
# FindSubCluster(, graph.name= "SCT_snn"), recognizes graph names with RNA_snn or SCT_snn in it.
# so, duplicate kmeans71 with a new name. 
# Algorithm for modularity optimization: 1 = original Louvain algorithm

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

DefaultAssay(sobj) <- 'SCT'

sobj$SCT_snn_km71 <- sobj$kmeans71

sobj <- SetIdent(sobj, value = "SCT_snn_km71")

# for the first cluster #1 
sobj <- FindSubCluster(sobj, cluster = "1", graph.name= "SCT_snn", 
                       subcluster.name = "sub71", resolution = 0.5, algorithm = 1)

# for the rest of the clusters
for (cluster in 2:71) {
  # Set the identity to the newly formed column 
  sobj <- SetIdent(sobj, value = "sub71")
  # Find subclusters 
  sobj <- FindSubCluster(sobj, cluster = as.character(cluster), 
                         graph.name = "SCT_snn", 
                         subcluster.name = "sub71", 
                         resolution = 0.5, algorithm = 1)
}

# tabulate and record sub-clusters by condition 

df <-as.data.frame.matrix(table(sobj$sub71, sobj$condition))

### calculation of normalized cell loss per sub-cluster. Only mDA sub-cluster were included in the final analysis. 

# Define beta value as FANS yield quotient (intact to lesioned)
beta <- 1.097  # see the article for explanation: Methods: Calculation of cell loss

# Calculate normalized_loss based on conditions per sub-cluster
df$normalized_loss <- ifelse(df$intact > df$lesion,
                             1 - ((df$lesion * beta) / df$intact),
                             (df$intact / (df$lesion * beta)) - 1)

df$subcluster <- rownames(df)

# write the name on the left of the "_" in subcluster as a new column cluster:
df$cluster <- gsub("\\_.*", "", df$subcluster)

# cell loss for only mDA neighborhoods and territories based on the sub-clusters. 
mda.clusters <- c('54', '65', '40', '13', '20', '52', '17', '28', '27', '45', '33', 
                  '23', '29', '66', '9', '67', '31', '44', '22', '46', '14', '11', 
                  '30', '1', '38', '41', '61', '21', '50', '19', '42', '4', '10', 
                  '2', '39', '60', '5', '43', '37', '53')

df <- df[df$cluster %in% mda.clusters, ]

all(unique(df$cluster) %in% mda.clusters )
#[1] TRUE

# set sub-clusters with negative cell loss to NA 
df$normalized_loss <- replace(df$normalized_loss, which(df$normalized_loss < 0), NA)

df <- df[, c(4, 5, 1:3)]
df <- na.omit(df)

write.csv(df, "/path/to/dir/cell_loss_km71subclusters.csv", row.names = TRUE)

# territory and neighborhood information (membership) was added to normalized cell loss per sub-clusters. 
# tests below were performed at cluster, neighborhood and territory levels;

# test of normality: 
shapiro.test(df$normalized_loss)
ggqqplot(df$normalized_loss, main= "only mDA sub-clusters loss")

# non-parametric Kruskal-Wallis rank sum test:
kruskal.test(df$normalized_loss ~ df$cluster, data = df)
kruskal.test(df$normalized_loss ~ df$Neighborhood, data = df)
kruskal.test(df$normalized_loss ~ df$Territory, data = df)

# post hoc pairwise tests: CI=the Conover-Iman test
CI.clusters <- as.data.frame(conover.test(df$normalized_loss, df$cluster, method="bh", list = TRUE))
CI.clusters <- CI.clusters[, c("comparisons", "P", "P.adjusted", "chi2", "T")]

CI.Neighborhood <- as.data.frame(conover.test(df$normalized_loss, df$Neighborhood, method="bh", list = TRUE))
CI.Neighborhood <- CI.Neighborhood[, c("comparisons", "P", "P.adjusted", "chi2", "T")]

CI.Territory <- as.data.frame(conover.test(df$normalized_loss, df$Territory, method="bh", list = TRUE))
CI.Territory <- CI.Territory[, c("comparisons", "P", "P.adjusted", "chi2", "T")]


### vulnerability / resilience modules genes selection:
# commonly up-regulated in intact-only nuclei of vulnerable / resilient clusters, 
# with p_val_adj < 0.05 & avg_log2FC > 0.5
# vulnerable clusters 54, 65, 13, 20, 45, 21 = normalized cell loss per cluster > 90%  
# resilient clusters  23, 29, 19, 42, 43, 53 = normalized cell loss per cluster < 50%  

# vulnerability & resilience modules 

vulmodule <- c('Kcnj6', 'Rgs6', 'Lrp1b', '9530059O14Rik', 'Lrrc3b', 'Nos1ap', 'Slc6a3', 'Cyyr1', 'Colgalt2', 'Ddc',
                 'Ripor2', 'Plcb4', 'Fam135b', 'Bnc2', 'Pou3f2', 'Pbx1', 'Klhl1', 'Gpr158', '1700042O10Rik', 'Vsnl1')

resmodule <- c('Kctd8', 'Frmd5', 'Atp8a1', 'Scg2', 'C130073E24Rik', 'Dtna', 'Lmx1a', 'Nav3')

# VUL module without DAT (Slc6a3)
NEWvulmodule <- c('Kcnj6', 'Rgs6', 'Lrp1b', '9530059O14Rik', 'Lrrc3b', 'Nos1ap', 'Cyyr1', 'Colgalt2', 'Ddc',
               'Ripor2', 'Plcb4', 'Fam135b', 'Bnc2', 'Pou3f2', 'Pbx1', 'Klhl1', 'Gpr158', '1700042O10Rik', 'Vsnl1')

# turn vectors into lists for adding as modules 

vulmodule <- list(vulmodule)

resmodule <- list(resmodule)

NEWvulmodule <- list(NEWvulmodule)

sobj <- AddModuleScore(object = sobj, features = vulmodule, name = "Vulmodule")

sobj <- AddModuleScore(object = sobj, features = resmodule, name = "Resmodule")

sobj <- AddModuleScore(object = sobj, features = NEWvulmodule, name = "NEWvulmodule")

# save the module scores as Assays.

sobj[['VUL']] <- CreateAssayObject(data = t(x = FetchData(object = sobj, vars = 'Vulmodule1')))

sobj[['RES']] <- CreateAssayObject(data = t(x = FetchData(object = sobj, vars = 'Resmodule1')))

sobj[['newVUL']] <- CreateAssayObject(data = t(x = FetchData(object = sobj, vars = 'NEWvulmodule1')))

# save the object 
saveRDS(sobj, "/path/to/dir/mDA.rds") 








sessionInfo()


