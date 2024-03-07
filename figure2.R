
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

### figure 2 + supplements 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

# Figure 2 UMAP 

p <- DimPlot(sobj, group.by = "kmeans71", cols = c("3"="#B35806", "56"="#B35806", "12"="#2171B5", "26"="#2171B5", "54"="#2171B5",
                                                   "65"="#2171B5", "40"="#2171B5", "13"="#2171B5", "20"="#2171B5", "52"="#2171B5", 
                                                   "17"="#2171B5", "28"="#2171B5", "27"="#2171B5", "45"="#2171B5", "33"="#2171B5", 
                                                   "23"="#2171B5", "29"="#2171B5", "66"="#2171B5", "9"="#2171B5", "67"="#2171B5", 
                                                   "31"="#2171B5", "44"="#2171B5", "22"="#2171B5", "46"="#2171B5", "14"="#2171B5", 
                                                   "11"="#2171B5","30"="#2171B5", "1"="#2171B5", "38"="#2171B5", "41"="#2171B5", 
                                                   "61"="#2171B5", "21"="#2171B5", "50"="#2171B5", "19"="#2171B5", "42"="#2171B5", 
                                                   "4"="#2171B5", "10"="#2171B5", "2"="#2171B5", "39"="#2171B5", "60"="#2171B5", 
                                                   "5"="#2171B5", "43"="#2171B5", "37"="#2171B5", "53"="#2171B5", "15"="#00441B", 
                                                   "71"="#00441B", "49"="#00441B", "68"="#00441B", "8"="#DF65B0", "62"="#DF65B0", 
                                                   "58"="#DF65B0", "47"="#DF65B0", "55"="#DF65B0", "16"="#DF65B0", "63"="#DF65B0", 
                                                   "6"="#DF65B0",  "7"="#DF65B0",  "24"="#DF65B0", "18"="#FB9A99", "57"="#FB9A99", 
                                                   "69"="#FB9A99", "48"="#FB9A99", "59"="#FB9A99", "70"="#FB9A99", "25"="#FB9A99", 
                                                   "34"="#FB9A99", "32"="#FB9A99", "64"="#FB9A99", "51"="#969696", "35"="#969696", 
                                                   "36"="#969696")) + coord_fixed() + NoLegend()

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=25, height=25, res=300)

LabelClusters(plot = p, id = 'kmeans71', color = 'white', size=4, fontface = 'bold', box = T, max.overlapp = Inf, repel=T) 

dev.off()


# Figure 2 dendrogram

dend <- as.dendrogram(hc71)

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=18, height=2, res=300)

p1 <- dend %>% set("leaves_pch", 15) %>%  
  set("leaves_cex", 0.01) %>%  
  set("branches_lwd", 0.4) %>%   
  set("leaves_col", c("3"="#B35806", "56"="#B35806", "12"="#2171B5", "26"="#2171B5", "54"="#2171B5",
                      "65"="#2171B5", "40"="#2171B5", "13"="#2171B5", "20"="#2171B5", "52"="#2171B5", 
                      "17"="#2171B5", "28"="#2171B5", "27"="#2171B5", "45"="#2171B5", "33"="#2171B5", 
                      "23"="#2171B5", "29"="#2171B5", "66"="#2171B5", "9"="#2171B5", "67"="#2171B5", 
                      "31"="#2171B5", "44"="#2171B5", "22"="#2171B5", "46"="#2171B5", "51"="#969696", 
                      "14"="#2171B5", "11"="#2171B5","30"="#2171B5", "1"="#2171B5", "38"="#2171B5", 
                      "41"="#2171B5", "61"="#2171B5", "21"="#2171B5", "50"="#2171B5", "19"="#2171B5", 
                      "42"="#2171B5", "4"="#2171B5", "10"="#2171B5", "2"="#2171B5", "39"="#2171B5", 
                      "60"="#2171B5", "5"="#2171B5", "43"="#2171B5", "37"="#2171B5", "53"="#2171B5", 
                      "8"="#DF65B0", "62"="#DF65B0", "15"="#00441B", "71"="#00441B", "49"="#00441B", 
                      "68"="#00441B", "47"="#DF65B0", "55"="#DF65B0", "16"="#DF65B0", "63"="#DF65B0", 
                      "36"="#969696", "35"="#969696", "6"="#DF65B0",  "7"="#DF65B0",  "24"="#DF65B0", 
                      "18"="#FB9A99", "57"="#FB9A99", "69"="#FB9A99", "48"="#FB9A99", "59"="#FB9A99", 
                      "58"="#DF65B0", "70"="#FB9A99", "25"="#FB9A99", "34"="#FB9A99", "32"="#FB9A99", 
                      "64"="#FB9A99"))

ggplot(p1, labels = F)

dev.off()


# Figure 2 dotplot

my.levels <- as.character(hc71$order)

# re-order levels according to dendrogram leaf nodes 

sobj$kmeans71  <- factor(x = sobj$kmeans71 , levels = my.levels)

markers <-  rev(unique(c("Mog", "Mag", 'Th', 'Slc6a3', 'Ddc', 'En1', 'Nr4a2', 'Pitx3', 
                         "Prlr", 'Satb2', "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Slc32a1")))  

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=55, height=20, res=300)

DotPlot(sobj, assay = "RNA", features = markers, group.by = "kmeans71") + coord_flip() + 
  theme(axis.text.x = element_text(size = 13, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(size = 18, face = "bold", angle = 0))

dev.off()

sessionInfo()
