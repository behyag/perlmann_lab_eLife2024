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

### figure 3 + supplements 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

#  Figure 3, Fig3_supplement: UMAP with neighborhoods & territories 

my.cols <- c("vector of customized colors")

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=25, height=25, res=300)
DimPlot(mda, group.by = "neighborhood", cols = my.cols, label = F) + coord_fixed() + 
  theme(legend.text=element_text(size = 12, face = 'bold'))

dev.off()

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=25, height=25, res=300)
DimPlot(mda, group.by = "territory", cols = my.cols, label = F) + coord_fixed() + 
  theme(legend.text=element_text(size = 12, face = 'bold'))

dev.off()

# Figure 3  DotPLot 

# the genes to be plotted for territories and neighborhoods: 

markers <- rev(unique(c("Th", "Slc6a3", "Zfp804b", "Mmp12", "Sprr1a", "Creb5", "Sox6", "Tigar", "Aldh1a1", "Vcan",
                        "Aldh1a7", "Anxa1", "Grin2c", "Wdr25", "Ndnf", "Calb1", "Slc17a6", "Slc32a1", "Gad2", "Ebf2", 
                        "Chrm2", "Zfp536", "Megf11", "Zeb2", "Met", "Fbn2", "Mkx", "Mid1", "Rxfp1", "Hs3st2", "Col23a1", 
                        "Dsg2", "Ism1", "Pcsk6", "Cald1", "Pde11a", "Tacr3", "Sema5b", "Pdia5", "Jph1", 
                        "Postn", "Arhgap28", "Npy1r", "Ebf1", "Cck", "Col24a1", "Npw", "Man1a", "Kctd8", "Vip", "Gipr", 
                        "Otx2", "Plpp4", "Plekhg1", "Grp", "Eya1", "Glra2", "Baiap3")))


# to put unassigned and non-mDA clusters on the far right of the plot, re-define the levels manually:

CL_levels <- c("12", "26", "54", "65", "40", "13", "20", "52", "17", "28", "27", "45", "33", "23", "29", "66",
               "9", "67", "31", "44", "22", "46", "51", "14", "11", "30", "1", "38", "41", "61", "21", "50", 
               "19", "42", "4", "10", "2", "39", "60", "5",  "43", "37", "53", "8", "62", "15", "71", "49", 
               "68", "47", "55", "16", "63", "36", "35", "6",  "7",  "24", "18", "57", "69", "48", "59", "58",
               "70", "25", "34", "32", "64", "3", "56")

sobj$kmeans71 <- factor(x = sobj$kmeans71, levels = CL_levels)

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=50, height=45, res=300)

DotPlot(sobj, assay = "RNA", features = markers, group.by = "kmeans71") + coord_flip() + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(size = 18, face = "bold", angle = 0))

dev.off()

# Sox6-Calb1 co-expression UMAP 

tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=75, height=25, res=300)

FeaturePlot(sobj, order = TRUE, cols = c("#EEEEEE", "blue", "orange"), pt.size = 0.5,
            features = c("Sox6", "Calb1"), blend = TRUE, blend.threshold = 0) + coord_fixed()

dev.off()

# Suppl Fig3  dotplot

complementary_markers <- rev(unique(c('Elavl2', 'Kcnj6', 'Nos1ap', 'Msi2', 'Pex5l', 'Lix1', 'Rgs6', 'Nwd2', 'Cdh11', 'Cntn4', 'Kcns3', 
                                      'Serpine2', 'Col25a1', 'Lmo3', 'Fam19a4', 'Col11a1', 'Vav3', 'Lama3', 'Vill', 'L3mbtl4', 
                                      'Zfp804b', 'Ntng1', 'Htr2c', 'Slc26a7', 'Crhbp', 'Pld1', 'Pde3a', 'Fgf10', 'Egfr', 'Slc17a6',    
                                      'Bend7', 'Bcl11a', 'Pard3b', 'Col14a1', 'Dkk2', 'Otx1', 'Pard3bos1', 'Gpc4', 'Lypd6b', 'Hpgd', 'Npnt',   
                                      'Trpc6', 'Igf1', 'Ranbp3l', 'Kank1', 'Cgnl1', '1700011I03Rik', 'Etv1', 'Atp8b1', 'Cpne2', 'Ano2',   
                                      'Adcy2', 'Bmpr1b', 'Cbln1', 'Htr2a', 'Lncenc1', 'Sulf1','Cpne5', 'Cyp26b1', 'Lpar1', 'Wnt7b',     
                                      'Mbnl3', 'Sorcs1', 'Fstl4', 'Syt9', 'Ppfibp2', 'Gulp1', 'Myo1b', 'Pou2f2',  
                                      'Kcnab1', 'Sorcs3', 'Otx2os1', 'Klhl14',  'Cbln4', 'Lpl', 'Neurod6', 'Grp', 'Piezo2', 'Fst', 
                                      'Nrp1', 'Pdgfd', 'Dab2', 'Adra1a', 'Nfib', 'Nfia', 'Pax5', 'Csf2rb2', 'Adcyap1')))


tiff(file = "/path/to/dir/plot_name.tiff", 
     units="cm", width=55, height=66, res=300)

DotPlot(sobj, assay = "RNA", features = complementary_markers, group.by = "kmeans71") + coord_flip() + 
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.text.y = element_text(size = 18, face = "bold", angle = 0))

dev.off()

sessionInfo()
