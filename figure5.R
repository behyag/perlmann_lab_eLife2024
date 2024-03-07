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
  library("graphics")
})

### figure 5 + supplements 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

sobj <- SetIdent(sobj, value = "condition")

tiff(file = "/path/to/dir/fig5_xx.tiff", 
     units="cm", width=40, height=20, res=300)

DimPlot(sobj, split.by = "condition", 
        cols = c("intact"="#08519C", "lesion"="#A63603")) + coord_fixed()

dev.off()



### Fig 5F: celltype dopaminergic vs non-dopaminergic

# class = cell type 

sobj <- SetIdent(sobj, value = "class")

sobj$class2 <- plyr::mapvalues(
  x = sobj$class, 
  from = c("midbrain_Dopaminergic", "GLUT", "GABA", "unassigned", "mODC", "Hy_DA"),
  to = rep(c('midbrain_Dopaminergic', 'non_Dopaminergic', 'Hy_DA'), c(1, 4, 1))) 

tb1 <- as.data.frame(table(sobj$class2, sobj$condition))

colnames(tb1)[1] <- 'celltype'
colnames(tb1)[2] <- 'condition'

# Stacked + percent

ggplot(tb1, aes(x=celltype, y=Freq , fill=condition )) +  
  geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(plot.title=element_text(size = 12, face = 'bold'), 
        legend.text = element_text(size = 10, face = 'bold'), 
        axis.text.x = element_text(size = 10, face = 'bold', angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10, face = 'bold')) +
  ggtitle('celltype by condition %') + 
  scale_y_continuous(expand = c(0,0)) 

# fig supplement 5 C Mosaic plot

# read in data with number of nuclei per condition per animal used 
d2 <- read_excel("/path/to/dir/condition_per_animal.xlsx", sheet = 1)

d2 <- as.data.frame(d2)

rownames(d2) <- d2$mouse_ID 

d2$mouse_ID <- NULL 

chsq <- chisq.test(d2)

mosaicplot(d2, shade = TRUE, las=2, cex.axis = 1.0, type = 'pearson',
           main = "condition per animal")


sessionInfo()
