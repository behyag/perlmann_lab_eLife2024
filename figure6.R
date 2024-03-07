suppressPackageStartupMessages({
  library(Seurat)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(future)
  require(scales)
  library(RColorBrewer)
  library("readxl")
  library(dplyr)
  library(dendextend)
  library('conover.test')
  library("devtools")
  library(ggpmisc)
})

### figure 6 + supplements 

# import the Dopaminergic nuclei (mDA) dataset.
sobj <- readRDS("/path/to/dir/mDA.rds")

# set levels by territory, but re-order again based on dendrogram (desired order) 
sobj <- SetIdent(sobj, value = "territory")

TER_levels <- c("ML", "Sox6", "Gad2", "Fbn2", ..., )

my.cols = c("vector of customized colors")

sobj$territory <- factor(x = sobj$territory, levels = TER_levels)

tiff(file = "/path/to/dir/plot.tiff", 
     units="cm", width=15, height=10, res=300)

VlnPlot(sobj, features = "Vulmodule1", cols = my.cols, group.by = "territory", pt.size = 0, combine = TRUE) +
  stat_summary(fun = mean, geom='point', size = 8, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "vulnerable module mDA territories")

dev.off()

tiff(file = "/path/to/dir/plot.tiff", 
     units="cm", width=15, height=10, res=300)

VlnPlot(sobj, features = "Resmodule1", cols = my.cols, group.by = "territory", pt.size = 0, combine = TRUE) +
  stat_summary(fun = mean, geom='point', size = 8, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "resilient module mDA territories")

dev.off()

# set levels by neighborhoods, but re-order again based on dendrogram (desired order) 
sobj <- SetIdent(sobj, value = "neighborhood")

NH_levels <- c("ML_NH1", "ML_NH2", "Sox6_NH1", "Sox6_NH2", ..., )

my.cols = c("vector of customized colors")

sobj$neighborhood <- factor(x = sobj$neighborhood, levels = NH_levels)

tiff(file = "/path/to/dir/plot.tiff", 
     units="cm", width=30, height=10, res=300)

VlnPlot(sobj, features = "Vulmodule1", cols = my.cols, group.by = "neighborhood", pt.size = 0, combine = TRUE) +
  stat_summary(fun = mean, geom='point', size = 8, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "vulnerable module mDA neighborhoods")

dev.off()

tiff(file = "/path/to/dir/plot.tiff", 
     units="cm", width=30, height=10, res=300)

VlnPlot(sobj, features = "Resmodule1", cols = my.cols, group.by = "neighborhood", pt.size = 0, combine = TRUE) +
  stat_summary(fun = mean, geom='point', size = 8, colour = "black", shape=95) + NoLegend() +
  labs(y="avg.Exp.gene.set", title = "resilient module mDA neighborhoods")

dev.off()


# Fig 6 E neighborhood pairwise comparison of Resilience module
# RES = resilience module 

UP.tb <- table(sobj@assays$RES@data, sobj$neighborhood)

UP.df <- as.data.frame(UP.tb)  

colnames(UP.df)
#[1] "Var1" "Var2" "Freq"

# neighborhoods re-ordered based on dendrogram and only choose mDA neighborhoods. 

mda.NH <- c('ML_NH1', 'ML_NH2', 'Sox6_NH1', 'Sox6_NH2', 'Sox6_NH3', 'Sox6_NH4', 
            'Gad2_NH1', 'Gad2_NH2', 'Fbn2_NH1', 'Fbn2_NH2', 'Pcsk6_NH1', 'Pcsk6_NH2', 
            'Pdia5_NH1', 'Pdia5_NH2', 'Col24a1', 'Vip', 'Otx2_NH1', 'Otx2_NH2')

# subset the defined neighborhoods in mda.NH from the UP.df 

for (i in mda.NH) {
  # Subset the dataframe for the current cluster
  cluster_df <- UP.df[UP.df$Var2 == i, ]
  
  # Filter the subsetted dataframe to retain only rows where Freq == 1
  cluster_df <- cluster_df[cluster_df$Freq == 1, ]
  
  # Assign the filtered dataframe to a new dataframe with the cluster name
  assign(i, cluster_df)
}

# merge into a long df:
ldf <- do.call('rbind', list(ML_NH1, ML_NH2, Sox6_NH1, Sox6_NH2, Sox6_NH3, Sox6_NH4, 
                             Gad2_NH1, Gad2_NH2, Fbn2_NH1, Fbn2_NH2, Pcsk6_NH1, Pcsk6_NH2, 
                             Pdia5_NH1, Pdia5_NH2, Col24a1, Vip, Otx2_NH1, Otx2_NH2))

colnames(ldf)
#[1] "Var1" "Var2" "Freq"

# large sample size, use lolcat package for test of normality
install_git("https://mikeburr.visualstudio.com/DefaultCollection/lolcat-public/_git/lolcat")

library(lolcat)

ldf$Var1 <- as.numeric(as.character(ldf$Var1 )) 

skewness.test(ldf$Var1)
# D'Agostino Skewness Normality Test

kurtosis.test(ldf$Var1 )
# D'Agostino Kurtosis Normality Test

ggqqplot(ldf$Var1 ) + ggtitle('vulnerable all NHs')

# test of homogeneity of variance
fligner.test(ldf$Var1 ~ ldf$Var2, data = ldf)

# non-parametric test (data is not normally distributed)

kruskal.test(ldf$Var1 ~ ldf$Var2, data = ldf)

# Welch's ANOVA because data shows heteroscedasticity (different groups have different standard deviations)

oneway.test(ldf$Var1 ~ ldf$Var2, data = ldf, var.equal = FALSE)

### conclusion: both tests reject the null hypothesis: there is significant difference among groups 

# post hoc pairwise tests: the Conover-Iman test
CI_NH_res <- as.data.frame(conover.test(ldf$Var1, ldf$Var2, method="bh", list = TRUE)) 

CI_NH_res <- CI_NH_res[, c("comparisons", "P", "P.adjusted", "chi2", "T")]

# visualization of pairwise comparisons for mDA neighborhoods: 
CI_NH_res$comparisons <- as.character(CI_NH_res$comparisons)

# remove white space in comparison column 
CI_NH_res$comparisons <- gsub('\\s+', '', CI_NH_res$comparisons)

# write the name on the right of the "-" in comparison string as a new column y:
CI_NH_res$y <- gsub(".*-", "", CI_NH_res$comparisons)

# write the name on the left of the "-" in comparison string as a new column x:
CI_NH_res$x <- gsub("\\-.*", "", CI_NH_res$comparisons)

# get the columns needed for visualization:

vis.df <- CI_NH_res[, c("x", "y", "P.adjusted")]

##  plot using geom_tile:
## for p <= alpha/2 (0.025)

tiff(file = "/path/to/dir/tilePlot.tiff", 
     units="cm", width=20, height=20, res=300)

ggplot(vis.df, aes(x, y, fill = P.adjusted)) + geom_tile() + 
  geom_hline(yintercept = seq_along(vis.df$y), color='grey') + 
  geom_vline(xintercept = seq_along(vis.df$x), color='grey') + 
  scale_fill_gradient(low = 'red', high = 'blue', limits=c(0, 0.025)) + 
  ggtitle('Res module comparison neighborhoods') + coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = 'bold')) +
  theme(axis.text.y = element_text(size = 14, face = 'bold')) 

dev.off()


# Fig Supplement 6 B  
# visualization of post hoc pairwise test: CI = Conover-Iman test
# CI.clusters derived from cell_loss.R script 

# remove white space in comparison column 
CI.clusters$comparisons <- gsub('\\s+', '', CI.clusters$comparisons)

# write the name on the right of the "-" in comparison string as a new column y:
CI.clusters$y <- gsub(".*-", "", CI.clusters$comparisons)

# write the name on the left of the "-" in comparison string as a new column x:
CI.clusters$x <- gsub("\\-.*", "", CI.clusters$comparisons)

# get the columns needed for visualization:

vis.df <- CI.clusters[, c("x", "y", "P.adjusted")]

## for p <= alpha/2 (0.025)
ggplot(vis.df, aes(x, y, fill = P.adjusted)) + geom_tile() + 
  geom_hline(yintercept = seq_along(vis.df$y), color='grey') + 
  geom_vline(xintercept = seq_along(vis.df$x), color='grey') + 
  scale_fill_gradient(low = 'red', high = 'blue', limits=c(0, 0.025)) + 
  ggtitle('normalized cell loss across mDA clusters')



### figure supplement 6 C & D 

# set cell loss as the dependent variable lm(Y~X) : Y: dependent   X:independent (predictor)
# set VUL (with DAT) or newVUL (without DAT) as predictor 
# calculate regression (linear model) fit linear regression models

# AverageExpression() Returns averaged expression values for each identity class
# Returns a matrix with genes as rows, identity classes as columns

vul <- as.data.frame(AverageExpression(
  object = sobj,
  assays = 'VUL',
  features = NULL,
  return.seurat = FALSE,
  group.by = "kmeans71",
  slot = "data",
  verbose = TRUE
))

newvul <- as.data.frame(AverageExpression(
  object = sobj,
  assays = 'newVUL',
  features = NULL,
  return.seurat = FALSE,
  group.by = "kmeans71",
  slot = "data",
  verbose = TRUE
))


## to drop the "vul" and "newVUL." from the column names

colnames(vul) <- gsub("VUL.", "", colnames(vul))

colnames(newvul) <- gsub("newVUL.", "", colnames(newvul))

# remove  non-mDA clusters, ML and unassigned clusters 12, 26, 51 
keepers <- c('40', '52', '17', '28', '27', '33', '23', '29', '66', '9', '67', '31', '44', 
              '22', '46', '14', '11', '30', '1', '38', '41', '61', '50', '19', '42', '4', '10', 
              '2', '39', '60', '5', '43', '37', '53')

vulmda <- vul[colnames(vul) %in% keepers] 

newvulmda <- newvul[colnames(newvul) %in% keepers] 

vulmda <- as.data.frame(t(vulmda))
newvulmda <- as.data.frame(t(newvulmda))

colnames(vulmda)[1] <- 'vul'
colnames(newvulmda)[1] <- 'newvul'

dfmda <- cbind.data.frame(vulmda, newvulmda)
dfmda$clusters <- rownames(dfmda)
dfmda <- dfmda[, c(3, 1, 2)]

### sqrt() moderate transformation to meet the normality assumption of linear models 

shapiro.test(sqrt(dfmda$vul))
# Shapiro-Wilk normality test
# data:  sqrt(dfmda$vul)
# W = 0.96427, p-value = 0.3224

shapiro.test(sqrt(dfmda$newvul))
# Shapiro-Wilk normality test
# data:  sqrt(dfmda$newvul)
# W = 0.96016, p-value = 0.2458

# from mDA.R script: load the normalized cell loss per mDA sub-clusters 
df <- read.csv(file = "path/to/file/cell_loss_km71subclusters.csv", header = TRUE, sep = ",", dec = ".")

# average cell loss per cluster
new_df <- df %>%
  group_by(cluster) %>%
  summarise(mean_normalized_loss = mean(normalized_loss))

# Get the order of clusters in dfmda (above)
order_clusters <- dfmda$clusters

# Reorder the rows of new_df to match the order of clusters in dfmda
new_df_reordered <- new_df[match(order_clusters, new_df$cluster), ]

new_df_reordered <- as.data.frame(new_df_reordered)

# create a new df for linear models 
lmdf <- cbind.data.frame(dfmda, new_df_reordered)

lmdf <- lmdf[, c(1:3, 5)]


### fit linear regression model for Vul (with DAT)

lm1 <- lm(mean_normalized_loss ~ sqrt(vul), data = lmdf)
summary(lm1)

# calculate a 95% confidence interval for the regression coefficient 

confint(lm1, 'sqrt(vul)', level = 0.95)
#                    2.5 %        97.5 %
#   sqrt(vul)    0.3468961     0.5536449


p1 <- ggplot(lm1, aes(mean_normalized_loss, lm1$model$`sqrt(vul)`)) + geom_point() +
  ggtitle("formula = avg.norm.loss ~ sqrt(Vul)") + stat_poly_eq() +
  geom_smooth(method="lm", col="red") + stat_regline_equation(label.x = 0, label.y = 1.25) +
  theme_classic()

p1 <- LabelPoints(plot = p1, points = colnames(lm1$residuals ),
                  size = 8,  color='blue', repel = T, xnudge = 0, ynudge = 0)
plot(p1)

### fit linear regression model for newVul (without DAT)

lm2 <- lm(mean_normalized_loss ~ sqrt(newvul), data = lmdf)
summary(lm2)

# calculate a 95% confidence interval for the regression coefficient

confint(lm2, 'sqrt(newvul)', level = 0.95)
#                    2.5 %     97.5 %
#   sqrt(newvul) 0.3592504  0.5927673


p2 <- ggplot(lm2, aes(mean_normalized_loss, lm2$model$`sqrt(newvul)`)) + geom_point() +
  ggtitle("formula = avg.norm.loss ~ sqrt(newVul)") + stat_poly_eq() +
  geom_smooth(method="lm", col="red") + stat_regline_equation(label.x = 0, label.y = 1.25) +
  theme_classic()

p2 <- LabelPoints(plot = p2, points = colnames(lm2$residuals ),
                  size = 8,  color='blue', repel = T, xnudge = 0, ynudge = 0)
plot(p2)



sessionInfo()


