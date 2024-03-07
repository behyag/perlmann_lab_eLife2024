
suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(sctransform)
  library(future)
  require(scales)
  library(RColorBrewer)
  library("readxl")
  library(dplyr)
  library(dendextend)
})


# After alignment with cellranger, (see cellrangercount.sh)
# the "filtered_feature_bc_matrix" per sample are read.

# Read sample info from CSV file
# the 12 sample IDs from P18856_3001 to P18856_3012. see supplementary file#1, sampleinfo_AllNuclei 
sample_info <- read.csv("/path/to/dir/sampleinfo_AllNuclei.csv")

# List of sample IDs
sample_ids <- sample_info$`10X_Serial_ID`

base_dir <- "/common/path/string/to/cellranger/output/files/"

# List to store Seurat objects
seurat_objs <- list()

# Iterate over each sample ID
for (sample_id in sample_ids) {
  # File path for filtered_feature_bc_matrix
  file_path <- paste0(base_dir, sample_id, "/outs/filtered_feature_bc_matrix")
  
  # Read in the data
  data <- Read10X(file_path)
  
  # Create Seurat objects
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_id)
  
  # Save individual Seurat object
  saveRDS(seurat_obj, file = paste0("/path/to/dir/", sample_id, ".rds"))
  
  # Add Seurat objects to the list
  seurat_objs[[sample_id]] <- seurat_obj
}

# Merge all Seurat objects into one object
merged_obj <- merge(x = seurat_objs[[1]], y = seurat_objs[-1], 
                    add.cell.ids = sample_ids)

# Calculate percentage of mitochondrial genes
merged_obj <- PercentageFeatureSet(merged_obj, "^mt-", col.name = "percent_mito")

# Calculate percentage of ribosomal genes
merged_obj <- PercentageFeatureSet(merged_obj, "^Rp[sl]", col.name = "percent_ribo")

# Cell cycle scoring
merged_obj <- CellCycleScoring(merged_obj, g2m.features = str_to_title(cc.genes$g2m.genes), 
                               s.features = str_to_title(cc.genes$s.genes))

# assign batch based on condition
assign_batch <- function(condition) {
  return(condition)
}

# Add batch to metadata using sample info
merged_obj$batch <- sapply(merged_obj$orig.ident, function(sample_id) {
  condition <- sample_info[sample_info$sample_ID == sample_id, "condition"]
  if (length(condition) > 0) {
    return(condition)
  } else {
    return("unknown")  # or any default value you prefer
  }
})

#  assign age based on age
assign_age <- function(age) {
  return(age)
}

# Add age to metadata using sample info
merged_obj$age <- sapply(merged_obj$orig.ident, function(sample_id) {
  age <- sample_info[sample_info$sample_ID == sample_id, "age"]
  if (length(age) > 0) {
    return(assign_age(age))
  } else {
    return("unknown")  # or any default value you prefer
  }
})

# Save merged Seurat object
saveRDS(merged_obj, file = "/path/to/dir/merged_seurat_obj.rds")


# (optional) load the individual saved objects for inspection & QC, etc.,  
# List of Seurat objects names
sample_ids <- c("s1", "s2", "s3", ...)

# Directory where Seurat objects are saved
dir_path <- "/path/to/dir/"

# List to store loaded Seurat objects
loaded_seurat_objects <- list()

# Iterate over each sample ID
for (sample_id in sample_ids) {
  # File path for the saved Seurat object
  file_path <- paste0(dir_path, sample_id, ".rds")
  
  # Load Seurat object
  seurat_obj <- readRDS(file_path)
  
  # Add Seurat object to the list
  loaded_seurat_objects[[sample_id]] <- seurat_obj
}

s1 <- loaded_seurat_objects[[1]]
s2 <- loaded_seurat_objects[[2]]
s3 <- loaded_seurat_objects[[3]]
#... etc.,


##  QC, pre-processing of the merged_obj and gene and cell filtering:

## cell filtering 
# remove cells with percent_mito > 5
# keep cells with nGenes within this range:  500 < nFeature_RNA < 10000 

## Gene filtering 
# keep genes which are expressed in 5 or more cells 

selected_cells <- WhichCells(merged_obj, expression = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent_mito < 5)

selected_genes <- rownames(merged_obj)[Matrix::rowSums(merged_obj) >= 10]

# the filtered Seurat object was defined as "sobj" for downstream analyses:
  
sobj <- subset(merged_obj, features = selected_genes, cells = selected_cells)

# remove Malat1

sobj <- sobj[ ! grepl("Malat1", rownames(sobj)), ]

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

# Cluster Determination:
sobj  <- FindClusters(sobj, resolution = seq(0.1, 0.9, by = 0.1), n.start = 100, n.iter = 100)

# rename clusters (1-26) instead of 0-25 for SCT_snn_res.0.1 resolution: 
sobj <- SetIdent(sobj, value = "SCT_snn_res.0.1")

newIDs <- rep(1:26, 1)

names(newIDs) <- levels(sobj)

sobj <- RenameIdents(sobj, newIDs) 

sobj@meta.data$SCT_snn_res.0.1 <- sobj@active.ident
  

# log-normalization genes in "RNA" assay 

# find 3000 HVGs for RNA assay so the same as in SCT assay. 

DefaultAssay(sobj) <- "RNA"

sobj <- NormalizeData(sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)


sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 3000)

# save the object 
saveRDS(sobj, "/path/to/dir/allnuc.rds") 


# cluster markers for Louvain resolution 0.1 

DefaultAssay(sobj) <- "RNA"

print("sobj")

sobj <- SetIdent(sobj, value = "SCT_snn_res.0.1")

print(length(levels(sobj@active.ident)))

plan("multicore", workers = 12)

options(future.globals.maxSize = 96000 * 1024^2)

for (i in 1:26){
  DE <- assign(paste0(i, "SCT_snn_res.0.1"), 
               FindMarkers(sobj, assay = "RNA", slot = "data", ident.1 = i, verbose = FALSE) )
  write.csv(DE, file=paste0("/path/to/dir/c_",i,".csv"))
}


sessionInfo()
