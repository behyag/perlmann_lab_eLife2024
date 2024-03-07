
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

# List of the 6 intact (6 mo) samples and 3 untreated young mice (3 mo) and 3 untreated old mice (18 mo). 
# for mice sample IDs see supplementary file#1, sampleinfo_intact_wt. 

# Read sample info from CSV file
sample_info <- read.csv("/path/to/dir/sampleinfo_intact_wt.csv")

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
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_id)
  
  # Add Seurat object to the list
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

# Define a function to assign batch based on condition
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



# Define a function to assign age based on age
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

### cell / gene filtering 

selected_cells <- WhichCells(merged_obj, expression = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent_mito < 5)

selected_genes <- rownames(merged_obj)[Matrix::rowSums(merged_obj) > 5]

sobj <- subset(merged_obj, features = selected_genes, cells = selected_cells)

# remove Malat1
sobj <- sobj[ ! grepl("Malat1", rownames(sobj)), ]

# dopaminergic nuclei filtering (mDA)
sobj <- subset(sobj, subset = Slc6a3 > 0 | Th > 0)

sobj <- SCTransform(sobj, variable.features.n = 1000)

sobj <- RunPCA(sobj, npcs = 100, verbose = FALSE)

ElbowPlot(sobj, reduction = "pca", ndims = 100)

sobj <- RunUMAP(sobj, dims = 1:30)

sobj <- RunTSNE(sobj, dims = 1:30)

DefaultAssay(sobj) <- "RNA"

sobj <- NormalizeData(sobj)

sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(sobj)

sobj <- ScaleData(sobj, features = all.genes)

#  PC 1:30   kmeans clustering  

pcmat <- Embeddings(sobj, reduction = 'pca')[,1:30]

set.seed(84)
km69s84 <- kmeans(pcmat, centers = 69, nstart = 50, iter.max = 1000, algorithm="MacQueen")

# add clusters to object
sobj@meta.data$km69s84 <- km69s84$cluster 

### identify and annotate clusters on this dataset based on markers' expression, dendrogram, etc, 

# create new metadata entry named "class" from annotated clusters 

sobj$class <- plyr::mapvalues(
  x = sobj$km69s84, 
  from = c('34', '9', '69', '6', '21', '13', '62', '39', '60', '7', '59', 
           '20', '63', '3', '12', '10', '19', '68', '1', '47', '31', '50', 
           '52', '23', '28', '53', '64', '65', '41', '2', '24', '27', '35', 
           '4', '42', '29', '44', '61', '55', '26', '37', '11', '8', '15', 
           '22', '32', '49', '14', '30', '51', '25', '46', '36', '67', '16', 
           '38', '56', '66', '18', '17', '40', '54', '58', '33', '5', '57', 
           '43', '45', '48'),
  to = rep(c('mDA', 'mODC', 'unassigned', 'Glut', 'Hy_DA', 'GABA', 'Glut'), 
           c(47, 1, 2, 4, 3, 10, 2)))

DimPlot(sobj, group.by = "class", cols = c("mODC" = "#B35806", "mDA" = "#2171B5", 
                                           "Glut" = "#DF65B0", "Hy_DA"="#00441B", 
                                           "GABA"="#FB9A99", "unassigned"="#969696"), 
        label = F, order = c('mODC', 'unassigned', 'GABA')) + coord_fixed() + 
  theme(text = element_text(size = 8, face = "bold"), legend.text=element_text(size = 14, face = 'bold'))

DimPlot(sobj, reduction = "umap", split.by = "batch", order = "untreated", 
        cols = c("intact"="#08519C", "untreated"="green")) + coord_fixed()

DefaultAssay(sobj) <- 'RNA'

markers <- c('Th', 'Slc6a3', 'Ddc', 'En1', 'Sox6', 'Calb1', 'Aldh1a1', 'Slc32a1', 'Gad2', 
             'Prlr', 'Satb2', 'Slc17a6', 'Nfib', 'Mog', 'Mag')

FeaturePlot(sobj, features = markers, slot = 'data', ncol = 5, coord.fixed = T)


sessionInfo()




