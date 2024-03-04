################################################################################
#                      Script for application on FIFA data                     #                                                                                 #
# In: Cavicchia, C. Vichi, M., Zaccaria, G. (2024+)                            #
# Parsimonious Gaussian mixture models. Revised version.                       #                                
################################################################################

# If necessary
# install.packages("datasetsICR")
library(datasetsICR)
# If necessary
# install.packages("mclust")
library(mclust)
# If necessary
# install.packages("vegan")
library(vegan)
# If necessary
# install.packages("dplyr")
library(dplyr)
# If necessary
# install.packages("tidyverse")
library(tidyverse)
# If necessary
# install.packages("remotes")
remotes::install_github("PUGMM-authors/PUGMM")
library(PUGMM)

# Loading the data
data(FIFA)
df <- FIFA
names <- df$Name
# PRE-PROCESSING
# Filtering out players with overall lower than or equal to 75
names <- names[df$Overall > 75]
df <- df[df$Overall > 75,]
df <- df[!duplicated(names),]
names <- names[!duplicated(names)]
# Filtering out goal keepers
names <- names[-which(df$Position == "GK")]
df <- df[-which(df$Position == "GK"),]
x <- df[, 46:74]

# ANALYSIS
# Running pugmm 
model <- pugmm(x, G = 2:10, m = 1:10)

# Displaying first 4 players per cluster
head(names[model$label == 1], 4)  # Cluster 1 :: strikers
head(names[model$label == 2], 4)  # Cluster 2 :: attacking midfielders
head(names[model$label == 3], 4)  # Cluster 3 :: stoppers
head(names[model$label == 4], 4)  # Cluster 4 :: central back
head(names[model$label == 5], 4)  # Cluster 5 :: number ten
head(names[model$label == 6], 4)  # Cluster 6 :: back wings
head(names[model$label == 7], 4)  # Cluster 7 :: midfielders
head(names[model$label == 8], 4)  # Cluster 8 :: forward

# Producing complete mds scatterplot
d <- dist(x)
cmd <- cmdscale(d)
groups <- levels(factor(model$label))
ordiplot(cmd, type = "n")
cols <- c("steelblue", "red", "green", "pink", "orange", "violet", "magenta", "grey")
for(i in seq_along(groups)){
  points(cmd[factor(model$label) == groups[i], ], col = cols[i], pch = 16)
}
ordispider(cmd, factor(model$label), label = TRUE)
ordihull(cmd, factor(model$label), lty = "dotted")

# Producing first 4 players per cluster mds scatterplot 
sample_size <- 4
cmds <- cbind(cmd, model$label, names)
xs <- matrix(rep(0, sample_size*8*4), sample_size*8, 4)
for (i in 1:8){
  xs[((i - 1)*sample_size + 1):(i*sample_size), ] <- as.matrix(head(cmds[model$label == i, ], sample_size))
}
groups <- levels(factor(model$label))
ordiplot(cmd, type = "n")
cols <- c("steelblue", "red", "green", "pink", "orange", "violet", "magenta", "grey")
for(i in seq_along(groups)){
  points(xs[factor(xs[, 3]) == groups[i], c(1:2)], col = cols[i], pch = 16, cex = .9)
}
text(xs[, c(1:2)], labels = xs[, 4], cex = .8)




