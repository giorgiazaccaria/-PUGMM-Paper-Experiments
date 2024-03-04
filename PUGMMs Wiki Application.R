################################################################################
#     Script for application on the use of Wikipedia within universities       #
#                                                                              #
# In: Cavicchia, C. Vichi, M., Zaccaria, G. (2024+)                            #
# Parsimonious Gaussian mixture models. Revised version.                       #                                
################################################################################

# If necessary
# install.packages("bnstruct")
library(bnstruct)
# If necessary
# install.packages("remotes")
remotes::install_github("PUGMM-authors/PUGMM")
library(PUGMM)

wiki <- read.csv("DIRECTORY  NAME/wiki.csv", sep = ";")

# PRE-PROCESSING
# Imputation of the missing variables
wiki.imp <- knn.impute(matrix(as.numeric(unlist(wiki)), dim(wiki)[1], dim(wiki)[2]), k = 5, cat.var = NULL)
colnames(wiki.imp) <- colnames(wiki)

# ANALYSIS
pugmm.wiki <- pugmm(wiki.imp, G = 1:10, m = 1:10)

pugmm.wiki$model.name
pugmm.wiki$G
pugmm.wiki$m
pugmm.wiki$V
pugmm.wiki$Sw
pugmm.wiki$Sb
plot.pugmm(pugmm.wiki)
