
suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(sctransform)
  library(future)
  require(scales)
  library(RColorBrewer)
  library(dplyr)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Matrix)
  library(ggplot2)
  library(harmony)
  library(WGCNA)
  library(igraph)
  library(devtools)
})

devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
packageVersion('hdWGCNA')
#  0.2.20

### Figure Supplement 6-3  co-expression network analysis 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

#  Set up Seurat object for WGCNA

sobj <- SetupForWGCNA(
  sobj,
  gene_select = "fraction", 
  fraction = 0.05, 
  wgcna_name = "con" # assign a name for the hdWGCNA experiment
)

# run harmony on PCA

pcmat <- Embeddings(sobj, reduction = 'pca')[,1:30]

harmonized_pcs <- HarmonyMatrix(
  data_mat  = pcmat,
  meta_data = sobj@meta.data ,
  vars_use  = "orig.ident",
  do_pca    = FALSE, 
  max.iter.harmony = 20,
  max.iter.cluster = 40,
  epsilon.harmony = -Inf, 
  epsilon.cluster = -Inf,
  plot_convergence = TRUE
)

## rename columns:
colnames(harmonized_pcs) <- paste0("H_", colnames(harmonized_pcs))

identical(rownames(harmonized_pcs), rownames(pcmat))
#  [1] TRUE

# Storing the custom dimensional reduction harmonized_pcs to Seurat reduction slot:
sobj[["harmony"]] <- CreateDimReducObject(embeddings = harmonized_pcs, key = "H_", assay = "RNA")

# construct metacells in each group but use 'harmony' instead of 'pca'
# construct metacells in each group on harmonized PCs

# groups are defined as group.by = 'territory.condition' 
sobj <- MetacellsByGroups(
  seurat_obj = sobj,
  group.by = 'territory.condition', 
  reduction = 'harmony', 
  k = 25, 
  max_shared = 10, 
  ident.group = 'territory.condition', 
  mode = "average",  
  min_cells = 100
)

###  Process the Metacell Seurat Object
# optional to get the metacell object from the hdWGCNA experiment using GetMetacellObject.
metacell_obj <- GetMetacellObject(sobj)

metacell_obj <- NormalizeData(metacell_obj)

metacell_obj <- FindVariableFeatures(metacell_obj)

metacell_obj <- ScaleData(metacell_obj, features = VariableFeatures(metacell_obj))

metacell_obj <- RunPCA(metacell_obj, features = VariableFeatures(metacell_obj))

ElbowPlot(metacell_obj, ndims = 50, reduction = "pca")  

sobj <- NormalizeMetacells(sobj)

sobj <- ScaleMetacells(sobj, features= VariableFeatures(sobj))

sobj <- RunPCAMetacells(sobj, features=VariableFeatures(sobj))

### 1:20 PCs were chosen based on metacell object elbowplot above
sobj <- RunUMAPMetacells(sobj, reduction='pca', dims=1:20)

# Figure Supplement 6-3 A 
DimPlotMetacells(sobj, group.by='territory3.condition', label=T, repel=T) + 
  umap_theme() + coord_fixed() + ggtitle("territory3.condition")


### Co-expression network analysis on the lesioned condition (mDA lesion)
## excluded ML-lesion and other non-mDA territories from network construction 
# Set up the expression matrix

sobj <- SetDatExpr(
  sobj,
  group_name = c('Ebf1_lesion', 'Fbn2_lesion', 'Gad2_lesion', 'Otx2_lesion',  
                 'Pcsk6_lesion', 'Pdia5_lesion', 'Sox6_lesion'),  
  group.by='territory.condition', # The same column used in MetacellsByGroups above 
  assay = 'RNA', 
  slot = 'data' 
)


# Select soft-power threshold
# Test different soft powers: Compute the scale-free topology model fit for different soft power thresholds

sobj <- TestSoftPowers(
  sobj,
  networkType = 'signed' 
)

power_table <- GetPowerTable(sobj)

write.csv(power_table, "/path/to/dir/soft_power_table.csv", row.names = F)

###  Construct co-expression network
sobj <- ConstructNetwork(
  sobj, soft_power=6,
  setDatExpr=FALSE,
  tom_name = 'lesionV2', 
  networkType = "signed",  
  TOMType = "signed", 
  overwrite_tom = TRUE
)

# visualize the WGCNA dendrogram, 
# grey module consists of genes that were not grouped into any co-expression module. 
# The grey module should be ignored for all downstream analysis and interpretations.

# Figure Supplement 6-3 B
PlotDendrogram(sobj, main='lesion V2 hdWGCNA Dendrogram')  

# hdWGCNA represents the co-expression network as a topoligcal overlap matrix (TOM). 
# This is a square matrix of genes by genes, where each value is the topoligcal overlap between the genes.  

TOM <- GetTOM(sobj)   # a gene x gene matrix

###  Module Eigengenes and Connectivity

sobj <- ScaleData(sobj, features=VariableFeatures(sobj))

sobj <- ModuleEigengenes(
  sobj,
  group.by.vars=NULL, 
  exclude_grey = TRUE, 
)

###  Compute module connectivity
sobj <- ModuleConnectivity(sobj, group.by = 'territory3.condition', 
                          group_name = c('Ebf1_lesion', 'Fbn2_lesion', 
                                         'Gad2_lesion', 'Otx2_lesion',
                                         'Pcsk6_lesion', 'Pdia5_lesion', 
                                         'Sox6_lesion'))

# rename the modules
sobj <- ResetModuleNames(
  sobj,
  new_name = "LmDAmod"
)

# optional to get the module assignment table:
modules <- GetModules(sobj)
write.csv(modules, "/path/to/dir/kMEs_df.csv", row.names = F)

# get the top 30 hub genes per module, sorted by kME 
hubdf <- GetHubGenes(sobj, n_hubs = 30)
write.csv(hubdf, "/path/to/dir/hub_df.csv", row.names = F)

#  visualize the correlation between each module based on their module eigengenes (MEs)
# Figure Supplement 6-3 E
ModuleCorrelogram(
  sobj,
  MEs2 = NULL,
  features = "MEs",
  order = "original",
  method = "ellipse",
  exclude_grey = TRUE,
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  sig.level =  0.05,
  pch.cex = 0.7,
  col = colorRampPalette(c("blue","white", "red"))(200),
  ncolors = 200,
  wgcna_name = NULL,
  wgcna_name2 = NULL
)

# plotting MEs 

# get MEs from seurat object
MEs <- GetMEs(sobj, harmonized=FALSE)

# add MEs to Seurat meta-data:
sobj@meta.data <- cbind(sobj@meta.data, MEs)

mods <- colnames(MEs); mods <- mods[mods != 'grey']

# re-order mods:

ord <- c("LmDAmod1", "LmDAmod2", "LmDAmod3", "LmDAmod4", 
         "LmDAmod5", "LmDAmod6", "LmDAmod7", "LmDAmod8", "LmDAmod9") 

mods <- mods[order(match(mods, ord))] 

p <- DotPlot(sobj, features=mods, group.by = 'territory')

# Figure Supplement 6-3 C
p + coord_flip() + RotatedAxis()

#### Network Visualization:  Combined hub gene network plots

# Figure Supplement 6-3 D

g <- HubGeneNetworkPlot(
  sobj,
  mods = 'all',
  n_hubs = 20,
  n_other = 5,
  sample_edges = TRUE,
  edge_prop = 0.65,
  return_graph = TRUE,
  edge.alpha = 0.95,
  vertex.label.cex = 0.2,
  hub.vertex.size = 3,
  other.vertex.size = 1,
  wgcna_name = NULL
)

# set seed before plotting with / without labels
set.seed(1357)
plot(g)

set.seed(1357)
plot(g, vertex.label=NA)

# Figure Supplement 6-3 F Boxplot kMEs of modules hub genes: 
# the colors for modules:
cols <- c('pink', 'blue', 'turquoise', 'yellow', 
          'brown', 'red', 'magenta', 'green', 'black')

ggplot(hubdf, aes(x=hubdf$module , y=hubdf$kME  )) + 
  geom_boxplot(width=0.45, fill=cols) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 14, face='bold', 
                                   vjust = 0.9, hjust= 0.9), 
        axis.text.y = element_text(size=14, face='bold'))


#### Gene Set Enrichment analysis

library(enrichR)
library(GeneOverlap)

theme_set(theme_cowplot())

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021',
         'GO_Cellular_Component_2021',
         'GO_Molecular_Function_2021')

# perform enrichment tests
sobj <- RunEnrichr(
  sobj,
  dbs=dbs, 
  max_genes = 100 
)

# retrieve the output table
enrich_df <- GetEnrichrTable(sobj)

# supplementary_file9 
write.csv(enrich_df, 
          "/path/to/dir/enrich_df.csv", 
          row.names = F)

### Figure Supplement 6-3 G enrichR plots for module 2

listEnrichrSites()

websiteLive <- TRUE

df <- enrich_df[ which( enrich_df$db == 'GO_Biological_Process_2021' & enrich_df$module == 'LmDAmod2') , ]

p <- if (websiteLive) plotEnrich(df, 
                                 showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value",
                                 title = 'GO_Biological_Process_2021_LmDAmodule2')

tiff(file = "/path/to/dir/module2_GO_BioProcess.tiff", 
     units="cm", width=22, height=10, res=300) 

p + theme(plot.title = element_text(size = 13, face = "bold"), 
          axis.text.x = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size = 10, face = "bold"))

dev.off()

df2 <- enrich_df[ which( enrich_df$db == 'GO_Biological_Process_2021' & enrich_df$module == 'LmDAmod3') , ]

p2 <- if (websiteLive) plotEnrich(df2, 
                                 showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value",
                                 title = 'GO_Biological_Process_2021_LmDAmodule3')

tiff(file = "/path/to/dir/module3_GO_BioProcess.tiff", 
     units="cm", width=22, height=10, res=300) 

p2 + theme(plot.title = element_text(size = 13, face = "bold"), 
          axis.text.x = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size = 10, face = "bold"))

dev.off()


sessionInfo()





