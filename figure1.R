
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

### figure 1 + figure supplements

# import all_nuclei dataset.
sobj <- readRDS("/path/to/dir/allnuc.rds")

# Seurat hierarchical clustering on the pseudobulk averages of different clusters

Idents(sobj) <- "SCT_snn_res.0.1"

sobj <- BuildClusterTree(sobj)

t1 <- Tool(sobj, slot = "BuildClusterTree")

ape::plot.phylo(t1, main = "AllCells hvg1000 SCT_snn_res.0.1", type = "phylogram", use.edge.length = TRUE,
                show.tip.label = TRUE, show.node.label = TRUE, edge.color = "black", edge.width = 1, edge.lty = 1,
                font = 3, cex = par("cex"), srt = 45, no.margin = FALSE, root.edge = TRUE, underscore = TRUE,
                direction = "downwards", tip.color = "black", plot = TRUE, align.tip.label=TRUE)


# turn tree (class=phylo) into a dendrogram object to use with dendextend package

dend <- as.dendrogram(t1)

# high resolution dendrogram dpi=300

tiff(file = "/path/to/dir/dend_name.tiff", 
     units="cm", width=5, height=25, res=300)

dend %>% set("leaves_pch", 15) %>%  
  set("leaves_cex", 3) %>%  
  set("labels_cex", 0.5) %>%   
  set("leaves_col", c("5"="#08519C", "23"="#253494", "25"="#00441B", "2"="#74C476", "11"="#C7E9C0",
                      "26"="#081D58", "16"="#B15928", "9"="#006D2C", "17"="#A1D99B", "6"="#6BAED6",
                      "20"="#9ECAE1", "13"="#E31A1C", "7" = "#E7D4E8", "10"="#C2A5CF", "1"="#762A83", 
                      "4"="#9970AB", "8"="#FD8D3C", "19"="#BD0026", "21"="#800026", "15"="#525252",  
                      "24"="#E7298A", "3"="#FDD0A2", "18"="#969696", "14"="#993404", 
                      "12"="#41AB5D", "22"="#FC4E2A")) %>% 
  plot(horiz=TRUE)  

dev.off()

# plot high res UMAP dpi=300

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=25, height=25, res=300)

DimPlot(ac, group.by = "LouvainRes.0.1", cols = c("1"="#762A83", "2"="#74C476", "3"="#FDD0A2",  
                                                  "4"="#9970AB", "5"="#08519C", "6"="#6BAED6", "7" = "#E7D4E8", 
                                                  "8"="#FD8D3C","9"="#006D2C", "10"="#C2A5CF", "11"="#C7E9C0", 
                                                  "12"="#41AB5D", "13"="#E31A1C", "14"="#993404", "15"="#525252", 
                                                  "16"="#B15928", "17"="#A1D99B", "18"="#969696", "19"="#BD0026", 
                                                  "20"="#9ECAE1", "21"="#800026", "22"="#FC4E2A", "23"="#253494", 
                                                  "24"="#E7298A", "25"="#00441B", "26"="#081D58"),     
        order = c("19", "24", "25"), label = T, label.size = 8, repel = T) + coord_fixed() + 
  theme(text = element_text(size = 8, face = "bold"), legend.text=element_text(size = 12, face = 'bold'))


dev.off()


# high res dotplot 300 dpi

# re-order object active levels (clusters) according to the tree leaf nodes. 

my.levels <- t1$tip.label  

sobj$SCT_snn_res.0.1  <- factor(x = sobj$SCT_snn_res.0.1 , levels = my.levels)

# import selected markers:
markers <- read_excel("/path/to/dir/Dot_plot_markers_Figxx.xlsx", col_names = T)

markers <- unique(pull(markers , Genes)) 

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=60, height=20, res=300)

DotPlot(sobj, assay = "RNA", features = markers, group.by = "LouvainRes.0.1") + 
  theme(axis.text.x = element_text(size = 15, face = "bold", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 15, face = "bold"))

dev.off()


sessionInfo()
