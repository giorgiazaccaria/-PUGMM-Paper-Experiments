################################################################################
#                      Script for benchmark data sets                          #
#                                                                              #
# In: Cavicchia, C. Vichi, M., Zaccaria, G. (2024+)                            #
# Parsimonious Gaussian mixture models. Revised version.                       #                                
################################################################################

# If necessary
# install.packages("MASS")
library(MASS)
# If necessary
# install.packages("mclust")
library(mclust)
# If necessary
# install.packages("pgmm")
library(pgmm)
# If necessary
# install.packages("HDclassif")
library(HDclassif)
# If necessary
# install.packages("IMIFA")
library(IMIFA)
# If necessary
# install.packages("remotes")
remotes::install_github("PUGMM-authors/PUGMM")
library(PUGMM)

################################ PENGUIN #######################################
penguin <- read.csv("DIRECTORY_NAME/penguins_size.txt", header = FALSE)
penguin <- penguin[, c(1, 3:6)]
penguin <- penguin[complete.cases(penguin), ]
x <- scale(penguin[, -1])
lab.penguin <- penguin[, 1]
G.max <- length(unique(lab.penguin)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.penguin <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.penguin <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.penguin <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.penguin <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.penguin$G
pugmm.penguin$m
pugmm.penguin$model.name
pugmm.penguin$pm.free
pugmm.penguin.ari <- adjustedRandIndex(lab.penguin, pugmm.penguin$label)
# GPCMs
mclust.penguin$G
mclust.penguin$modelName
nMclustParams(mclust.penguin$modelName, mclust.penguin$d, mclust.penguin$G)
mclust.penguin.ari <- adjustedRandIndex(lab.penguin, mclust.penguin$classification)
# PGMMs
pgmm.penguin$g
pgmm.penguin$q
pgmm.penguin$model
PGMM_dfree(pgmm.penguin$q, mclust.penguin$d, pgmm.penguin$g, pgmm.penguin$model)
pgmm.penguin.ari <- adjustedRandIndex(lab.penguin, pgmm.penguin$map)
# HDDC
hddc.penguin$K
hddc.penguin$d
hddc.penguin$model
hddc.penguin$complexity
hddc.penguin.ari <- adjustedRandIndex(lab.penguin, hddc.penguin$class)

# Gopt
# PUGMMs
if(pugmm.penguin$G != length(unique(lab.penguin))){
  pugmm.penguin.Gopt <- pugmm(x, G = length(unique(lab.penguin)), m = 1:m.max)
  pugmm.penguin.Gopt$m
  pugmm.penguin.Gopt$model.name
  pugmm.penguin.Gopt$pm.free
  pugmm.penguin.Gopt.ari <- adjustedRandIndex(lab.penguin, pugmm.penguin.Gopt$label)}
# GPCMs
if(mclust.penguin$G != length(unique(lab.penguin))){
  mclust.penguin.Gopt <- Mclust(x, G = length(unique(lab.penguin)), control = emControl(itmax = 500))
  mclust.penguin.Gopt$modelName
  nMclustParams(mclust.penguin.Gopt$modelName, mclust.penguin.Gopt$d, mclust.penguin.Gopt$G)
  mclust.penguin.Gopt.ari <- adjustedRandIndex(lab.penguin, mclust.penguin.Gopt$classification)}
# PGMMs
if(pgmm.penguin$g != length(unique(lab.penguin))){
  pgmm.penguin.Gopt <- pgmmEM(x, rG = length(unique(lab.penguin)), rq = 1:m.max, relax = TRUE)
  pgmm.penguin.Gopt$q
  pgmm.penguin.Gopt$model
  PGMM_dfree(pgmm.penguin.Gopt$q, mclust.penguin$d, pgmm.penguin.Gopt$g, pgmm.penguin.Gopt$model) 
  pgmm.penguin.Gopt.ari <- adjustedRandIndex(lab.penguin, pgmm.penguin.Gopt$map)}
# HDDC
if(hddc.penguin$K != length(unique(lab.penguin))){
  hddc.penguin.Gopt <- hddc(x, K = length(unique(lab.penguin)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.penguin.Gopt$model
  hddc.penguin.Gopt$complexity
  hddc.penguin.Gopt.ari <- adjustedRandIndex(lab.penguin, hddc.penguin.Gopt$class)}

################################ WINE 13 #######################################
data(wine, package = "HDclassif")
x <- scale(wine[, -1])
lab.wine <- wine[, 1]
G.max <- length(unique(lab.wine)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.wine <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.wine <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.wine <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.wine <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.wine$G
pugmm.wine$m
pugmm.wine$model.name
pugmm.wine$pm.free
pugmm.wine.ari <- adjustedRandIndex(lab.wine, pugmm.wine$label)
# GPCMs
mclust.wine$G
mclust.wine$modelName
nMclustParams(mclust.wine$modelName, mclust.wine$d, mclust.wine$G)
mclust.wine.ari <- adjustedRandIndex(lab.wine, mclust.wine$classification)
# PGMMs
pgmm.wine$g
pgmm.wine$q
pgmm.wine$model
PGMM_dfree(pgmm.wine$q, mclust.wine$d, pgmm.wine$g, pgmm.wine$model)
pgmm.wine.ari <- adjustedRandIndex(lab.wine, pgmm.wine$map)
# HDDC
hddc.wine$K
hddc.wine$d
hddc.wine$model
hddc.wine$complexity
hddc.wine.ari <- adjustedRandIndex(lab.wine, hddc.wine$class)

# Gopt
# PUGMMs
if(pugmm.wine$G != length(unique(lab.wine))){
  pugmm.wine.Gopt <- pugmm(x, G = length(unique(lab.wine)), m = 1:m.max)
  pugmm.wine.Gopt$m
  pugmm.wine.Gopt$model.name
  pugmm.wine.Gopt$pm.free
  pugmm.wine.Gopt.ari <- adjustedRandIndex(lab.wine, pugmm.wine.Gopt$label)}
# GPCMs
if(mclust.wine$G != length(unique(lab.wine))){
  mclust.wine.Gopt <- Mclust(x, G = length(unique(lab.wine)), control = emControl(itmax = 500))
  mclust.wine.Gopt$modelName
  nMclustParams(mclust.wine.Gopt$modelName, mclust.wine.Gopt$d, mclust.wine.Gopt$G)
  mclust.wine.Gopt.ari <- adjustedRandIndex(lab.wine, mclust.wine.Gopt$classification)}
# PGMMs
if(pgmm.wine$g != length(unique(lab.wine))){
  pgmm.wine.Gopt <- pgmmEM(x, rG = length(unique(lab.wine)), rq = 1:m.max, relax = TRUE)
  pgmm.wine.Gopt$q
  pgmm.wine.Gopt$model
  PGMM_dfree(pgmm.wine.Gopt$q, mclust.wine$d, pgmm.wine.Gopt$g, pgmm.wine.Gopt$model) 
  pgmm.wine.Gopt.ari <- adjustedRandIndex(lab.wine, pgmm.wine.Gopt$map)}
# HDDC
if(hddc.wine$K != length(unique(lab.wine))){
  hddc.wine.Gopt <- hddc(x, K = length(unique(lab.wine)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.wine.Gopt$model
  hddc.wine.Gopt$complexity
  hddc.wine.Gopt.ari <- adjustedRandIndex(lab.wine, hddc.wine.Gopt$class)}

################################ WINE 27 #######################################
data(wine, package = "pgmm")
x <- scale(wine[, -1])
lab.wine27 <- wine[, 1]
G.max <- length(unique(lab.wine27)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.wine27 <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.wine27 <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.wine27 <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.wine27 <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.wine27$G
pugmm.wine27$m
pugmm.wine27$model.name
pugmm.wine27$pm.free
pugmm.wine27.ari <- adjustedRandIndex(lab.wine27, pugmm.wine27$label)
# GPCMs
mclust.wine27$G
mclust.wine27$modelName
nMclustParams(mclust.wine27$modelName, mclust.wine27$d, mclust.wine27$G)
mclust.wine27.ari <- adjustedRandIndex(lab.wine27, mclust.wine27$classification)
# PGMMs
pgmm.wine27$g
pgmm.wine27$q
pgmm.wine27$model
PGMM_dfree(pgmm.wine27$q, mclust.wine27$d, pgmm.wine27$g, pgmm.wine27$model)
pgmm.wine27.ari <- adjustedRandIndex(lab.wine27, pgmm.wine27$map)
# HDDC
hddc.wine27$K
hddc.wine27$d
hddc.wine27$model
hddc.wine27$complexity
hddc.wine27.ari <- adjustedRandIndex(lab.wine27, hddc.wine27$class)

# Gopt
# PUGMMs
if(pugmm.wine27$G != length(unique(lab.wine27))){
  pugmm.wine27.Gopt <- pugmm(x, G = length(unique(lab.wine27)), m = 1:m.max)
  pugmm.wine27.Gopt$m
  pugmm.wine27.Gopt$model.name
  pugmm.wine27.Gopt$pm.free
  pugmm.wine27.Gopt.ari <- adjustedRandIndex(lab.wine27, pugmm.wine27.Gopt$label)}
# GPCMs
if(mclust.wine27$G != length(unique(lab.wine27))){
  mclust.wine27.Gopt <- Mclust(x, G = length(unique(lab.wine27)), control = emControl(itmax = 500))
  mclust.wine27.Gopt$modelName
  nMclustParams(mclust.wine27.Gopt$modelName, mclust.wine27.Gopt$d, mclust.wine27.Gopt$G)
  mclust.wine27.Gopt.ari <- adjustedRandIndex(lab.wine27, mclust.wine27.Gopt$classification)}
# PGMMs
if(pgmm.wine27$g != length(unique(lab.wine27))){
  pgmm.wine27.Gopt <- pgmmEM(x, rG = length(unique(lab.wine27)), rq = 1:m.max, relax = TRUE)
  pgmm.wine27.Gopt$q
  pgmm.wine27.Gopt$model
  PGMM_dfree(pgmm.wine27.Gopt$q, mclust.wine27$d, pgmm.wine27.Gopt$g, pgmm.wine27.Gopt$model) 
  pgmm.wine27.Gopt.ari <- adjustedRandIndex(lab.wine27, pgmm.wine27.Gopt$map)}
# HDDC
if(hddc.wine27$K != length(unique(lab.wine27))){
  hddc.wine27.Gopt <- hddc(x, K = length(unique(lab.wine27)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.wine27.Gopt$model
  hddc.wine27.Gopt$complexity
  hddc.wine27.Gopt.ari <- adjustedRandIndex(lab.wine27, hddc.wine27.Gopt$class)}

################################ THYROID #######################################
data(thyroid, package = "mclust")
thyroid <- as.data.frame(thyroid)
lab.thyroid <- as.numeric(thyroid[, 1])
x <- scale(thyroid[, -1])
G.max <- length(unique(lab.thyroid)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.thyroid <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.thyroid <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.thyroid <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.thyroid <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.thyroid$G
pugmm.thyroid$m
pugmm.thyroid$model.name
pugmm.thyroid$pm.free
pugmm.thyroid.ari <- adjustedRandIndex(lab.thyroid, pugmm.thyroid$label)
# GPCMs
mclust.thyroid$G
mclust.thyroid$modelName
nMclustParams(mclust.thyroid$modelName, mclust.thyroid$d, mclust.thyroid$G)
mclust.thyroid.ari <- adjustedRandIndex(lab.thyroid, mclust.thyroid$classification)
# PGMMs
# pgmm.thyroid$g
# pgmm.thyroid$q
# pgmm.thyroid$model                                                                   
# PGMM_dfree(pgmm.thyroid$q, mclust.thyroid$d, pgmm.thyroid$g, pgmm.thyroid$model)
# pgmm.thyroid.ari <- adjustedRandIndex(lab.thyroid, pgmm.thyroid$map)
# HDDC
hddc.thyroid$K
hddc.thyroid$d
hddc.thyroid$model
hddc.thyroid$complexity
hddc.thyroid.ari <- adjustedRandIndex(lab.thyroid, hddc.thyroid$class)

# Gopt
# PUGMMs
if(pugmm.thyroid$G != length(unique(lab.thyroid))){
  pugmm.thyroid.Gopt <- pugmm(x, G = length(unique(lab.thyroid)), m = 1:m.max)
  pugmm.thyroid.Gopt$m
  pugmm.thyroid.Gopt$model.name
  pugmm.thyroid.Gopt$pm.free
  pugmm.thyroid.Gopt.ari <- adjustedRandIndex(lab.thyroid, pugmm.thyroid.Gopt$label)}
# GPCMs
if(mclust.thyroid$G != length(unique(lab.thyroid))){
  mclust.thyroid.Gopt <- Mclust(x, G = length(unique(lab.thyroid)), control = emControl(itmax = 500))
  mclust.thyroid.Gopt$modelName
  nMclustParams(mclust.thyroid.Gopt$modelName, mclust.thyroid.Gopt$d, mclust.thyroid.Gopt$G)
  mclust.thyroid.Gopt.ari <- adjustedRandIndex(lab.thyroid, mclust.thyroid.Gopt$classification)}
# PGMMs
pgmm.thyroid.Gopt <- pgmmEM(x, rG = length(unique(lab.thyroid)), rq = 1:m.max, relax = TRUE)
pgmm.thyroid.Gopt$q
pgmm.thyroid.Gopt$model
PGMM_dfree(pgmm.thyroid.Gopt$q, mclust.thyroid$d, pgmm.thyroid.Gopt$g, pgmm.thyroid.Gopt$model) 
pgmm.thyroid.Gopt.ari <- adjustedRandIndex(lab.thyroid, pgmm.thyroid.Gopt$map)
# HDDC
if(hddc.thyroid$K != length(unique(lab.thyroid))){
  hddc.thyroid.Gopt <- hddc(x, K = length(unique(lab.thyroid)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.thyroid.Gopt$model
  hddc.thyroid.Gopt$complexity
  hddc.thyroid.Gopt.ari <- adjustedRandIndex(lab.thyroid, hddc.thyroid.Gopt$class)}

################################################################################
################################ KIDNEY ########################################
data(ckd, package = "teigen")
x <- scale(ckd[, -1])
lab.kidney <- as.numeric(ckd[, 1])
G.max <- length(unique(lab.kidney)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.kidney <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.kidney <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.kidney <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.kidney <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.kidney$G
pugmm.kidney$m
pugmm.kidney$model.name
pugmm.kidney$pm.free
pugmm.kidney.ari <- adjustedRandIndex(lab.kidney, pugmm.kidney$label)
# GPCMs
mclust.kidney$G
mclust.kidney$modelName
nMclustParams(mclust.kidney$modelName, mclust.kidney$d, mclust.kidney$G)
mclust.kidney.ari <- adjustedRandIndex(lab.kidney, mclust.kidney$classification)
# PGMMs
pgmm.kidney$g
pgmm.kidney$q
pgmm.kidney$model
PGMM_dfree(pgmm.kidney$q, mclust.kidney$d, pgmm.kidney$g, pgmm.kidney$model)
pgmm.kidney.ari <- adjustedRandIndex(lab.kidney, pgmm.kidney$map)
# HDDC
hddc.kidney$K
hddc.kidney$d
hddc.kidney$model
hddc.kidney$complexity
hddc.kidney.ari <- adjustedRandIndex(lab.kidney, hddc.kidney$class)

# Gopt
# PUGMMs
if(pugmm.kidney$G != length(unique(lab.kidney))){
  pugmm.kidney.Gopt <- pugmm(x, G = length(unique(lab.kidney)), m = 1:m.max)
  pugmm.kidney.Gopt$m
  pugmm.kidney.Gopt$model.name
  pugmm.kidney.Gopt$pm.free
  pugmm.kidney.Gopt.ari <- adjustedRandIndex(lab.kidney, pugmm.kidney.Gopt$label)}
# GPCMs
if(mclust.kidney$G != length(unique(lab.kidney))){
  mclust.kidney.Gopt <- Mclust(x, G = length(unique(lab.kidney)), control = emControl(itmax = 500))
  mclust.kidney.Gopt$modelName
  nMclustParams(mclust.kidney.Gopt$modelName, mclust.kidney.Gopt$d, mclust.kidney.Gopt$G)
  mclust.kidney.Gopt.ari <- adjustedRandIndex(lab.kidney, mclust.kidney.Gopt$classification)}
# PGMMs
if(pgmm.kidney$g != length(unique(lab.kidney))){
  pgmm.kidney.Gopt <- pgmmEM(x, rG = length(unique(lab.kidney)), rq = 1:m.max, relax = TRUE)
  pgmm.kidney.Gopt$q
  pgmm.kidney.Gopt$model
  PGMM_dfree(pgmm.kidney.Gopt$q, mclust.kidney$d, pgmm.kidney.Gopt$g, pgmm.kidney.Gopt$model)
  pgmm.kidney.Gopt.ari <- adjustedRandIndex(lab.kidney, pgmm.kidney.Gopt$map)}
# HDDC
if(hddc.kidney$K != length(unique(lab.kidney))){
  hddc.kidney.Gopt <- hddc(x, K = length(unique(lab.kidney)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.kidney.Gopt$model
  hddc.kidney.Gopt$complexity
  hddc.kidney.Gopt.ari <- adjustedRandIndex(lab.kidney, hddc.kidney.Gopt$class)}

################################ ECONOMICS #####################################
data(Economics, package = "datasetsICR")
x <- scale(Economics[, -13])
lab.economics <- as.numeric(Economics[, 13])
G.max <- length(unique(lab.economics)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.economics <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.economics <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.economics <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.economics <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.economics$G
pugmm.economics$m
pugmm.economics$model.name
pugmm.economics$pm.free
pugmm.economics.ari <- adjustedRandIndex(lab.economics, pugmm.economics$label)
# GPCMs
mclust.economics$G
mclust.economics$modelName
nMclustParams(mclust.economics$modelName, mclust.economics$d, mclust.economics$G)
mclust.economics.ari <- adjustedRandIndex(lab.economics, mclust.economics$classification)
# PGMMs
# pgmm.economics$g
# pgmm.economics$q
# pgmm.economics$model                                                                   ERROR OCCURED               
# PGMM_dfree(pgmm.economics$q, mclust.economics$d, pgmm.economics$g, pgmm.economics$model)
# pgmm.economics.ari <- adjustedRandIndex(lab.economics, pgmm.economics$map)
# HDDC
hddc.economics$K
hddc.economics$d
hddc.economics$model
hddc.economics$complexity
hddc.economics.ari <- adjustedRandIndex(lab.economics, hddc.economics$class)

# Gopt
# PUGMMs
if(pugmm.economics$G != length(unique(lab.economics))){
  pugmm.economics.Gopt <- pugmm(x, G = length(unique(lab.economics)), m = 1:m.max)
  pugmm.economics.Gopt$m
  pugmm.economics.Gopt$model.name
  pugmm.economics.Gopt$pm.free
  pugmm.economics.Gopt.ari <- adjustedRandIndex(lab.economics, pugmm.economics.Gopt$label)}
# GPCMs
if(mclust.economics$G != length(unique(lab.economics))){
  mclust.economics.Gopt <- Mclust(x, G = length(unique(lab.economics)), control = emControl(itmax = 500))
  mclust.economics.Gopt$modelName
  nMclustParams(mclust.economics.Gopt$modelName, mclust.economics.Gopt$d, mclust.economics.Gopt$G)
  mclust.economics.Gopt.ari <- adjustedRandIndex(lab.economics, mclust.economics.Gopt$classification)}
# PGMMs
pgmm.economics.Gopt <- pgmmEM(x, rG = length(unique(lab.economics)), rq = 1:m.max, relax = TRUE)
pgmm.economics.Gopt$q
pgmm.economics.Gopt$model
PGMM_dfree(pgmm.economics.Gopt$q, mclust.economics$d, pgmm.economics.Gopt$g, pgmm.economics.Gopt$model) 
pgmm.economics.Gopt.ari <- adjustedRandIndex(lab.economics, pgmm.economics.Gopt$map)
# HDDC
if(hddc.economics$K != length(unique(lab.economics))){
  hddc.economics.Gopt <- hddc(x, K = length(unique(lab.economics)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.economics.Gopt$model
  hddc.economics.Gopt$complexity
  hddc.economics.Gopt.ari <- adjustedRandIndex(lab.economics, hddc.economics.Gopt$class)}

################################ TETRAGONULA ###################################
tetragonula <- read.csv("DIRECTORY_NAME/Tetragonula.txt", header = TRUE, sep = ";", dec = ",")
x <- scale(apply(tetragonula[, 3:6], 2, as.numeric))
lab.tetragonula <- tetragonula[, 2]
G.max <- length(unique(lab.tetragonula)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.tetragonula <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.tetragonula <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.tetragonula <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.tetragonula <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.tetragonula$G
pugmm.tetragonula$m
pugmm.tetragonula$model.name
pugmm.tetragonula$pm.free
pugmm.tetragonula.ari <- adjustedRandIndex(lab.tetragonula, pugmm.tetragonula$label)
# GPCMs
mclust.tetragonula$G
mclust.tetragonula$modelName
nMclustParams(mclust.tetragonula$modelName, mclust.tetragonula$d, mclust.tetragonula$G)
mclust.tetragonula.ari <- adjustedRandIndex(lab.tetragonula, mclust.tetragonula$classification)
# PGMMs
pgmm.tetragonula$g
pgmm.tetragonula$q
pgmm.tetragonula$model                                                                  
PGMM_dfree(pgmm.tetragonula$q, mclust.tetragonula$d, pgmm.tetragonula$g, pgmm.tetragonula$model)
pgmm.tetragonula.ari <- adjustedRandIndex(lab.tetragonula, pgmm.tetragonula$map)
# HDDC
hddc.tetragonula$K
hddc.tetragonula$d
hddc.tetragonula$model
hddc.tetragonula$complexity
hddc.tetragonula.ari <- adjustedRandIndex(lab.tetragonula, hddc.tetragonula$class)

# Gopt
# PUGMMs
if(pugmm.tetragonula$G != length(unique(lab.tetragonula))){
  pugmm.tetragonula.Gopt <- pugmm(x, G = length(unique(lab.tetragonula)), m = 1:m.max)
  pugmm.tetragonula.Gopt$m
  pugmm.tetragonula.Gopt$model.name
  pugmm.tetragonula.Gopt$pm.free
  pugmm.tetragonula.Gopt.ari <- adjustedRandIndex(lab.tetragonula, pugmm.tetragonula.Gopt$label)}
# GPCMs
if(mclust.tetragonula$G != length(unique(lab.tetragonula))){
  mclust.tetragonula.Gopt <- Mclust(x, G = length(unique(lab.tetragonula)), control = emControl(itmax = 500))
  mclust.tetragonula.Gopt$modelName
  nMclustParams(mclust.tetragonula.Gopt$modelName, mclust.tetragonula.Gopt$d, mclust.tetragonula.Gopt$G)
  mclust.tetragonula.Gopt.ari <- adjustedRandIndex(lab.tetragonula, mclust.tetragonula.Gopt$classification)}
# PGMMs
if(pgmm.tetragonula$g != length(unique(lab.tetragonula))){
  pgmm.tetragonula.Gopt <- pgmmEM(x, rG = length(unique(lab.tetragonula)), rq = 1:m.max, relax = TRUE)
  pgmm.tetragonula.Gopt$q
  pgmm.tetragonula.Gopt$model
  PGMM_dfree(pgmm.tetragonula.Gopt$q, mclust.tetragonula$d, pgmm.tetragonula.Gopt$g, pgmm.tetragonula.Gopt$model) 
  pgmm.tetragonula.Gopt.ari <- adjustedRandIndex(lab.tetragonula, pgmm.tetragonula.Gopt$map)}
# HDDC
if(hddc.tetragonula$K != length(unique(lab.tetragonula))){
  hddc.tetragonula.Gopt <- hddc(x, K = length(unique(lab.tetragonula)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.tetragonula.Gopt$model
  hddc.tetragonula.Gopt$complexity
  hddc.tetragonula.Gopt.ari <- adjustedRandIndex(lab.tetragonula, hddc.tetragonula.Gopt$class)}

################################ DIABETES ######################################
data(diabetes, package = "mclust")
x <- scale(diabetes[, -1])
lab.diabetes <- as.numeric(diabetes[, 1])
G.max <- length(unique(lab.diabetes)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.diabetes <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.diabetes <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.diabetes <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.diabetes <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.diabetes$G
pugmm.diabetes$m
pugmm.diabetes$model.name
pugmm.diabetes$pm.free
pugmm.diabetes.ari <- adjustedRandIndex(lab.diabetes, pugmm.diabetes$label)
# GPCMs
mclust.diabetes$G
mclust.diabetes$modelName
nMclustParams(mclust.diabetes$modelName, mclust.diabetes$d, mclust.diabetes$G)
mclust.diabetes.ari <- adjustedRandIndex(lab.diabetes, mclust.diabetes$classification)
# PGMMs
pgmm.diabetes$g
pgmm.diabetes$q
pgmm.diabetes$model
PGMM_dfree(pgmm.diabetes$q, mclust.diabetes$d, pgmm.diabetes$g, pgmm.diabetes$model)
pgmm.diabetes.ari <- adjustedRandIndex(lab.diabetes, pgmm.diabetes$map)
# HDDC
hddc.diabetes$K
hddc.diabetes$d
hddc.diabetes$model
hddc.diabetes$complexity
hddc.diabetes.ari <- adjustedRandIndex(lab.diabetes, hddc.diabetes$class)

# Gopt
# PUGMMs
if(pugmm.diabetes$G != length(unique(lab.diabetes))){
  pugmm.diabetes.Gopt <- pugmm(x, G = length(unique(lab.diabetes)), m = 1:m.max)
  pugmm.diabetes.Gopt$m
  pugmm.diabetes.Gopt$model.name
  pugmm.diabetes.Gopt$pm.free
  pugmm.diabetes.Gopt.ari <- adjustedRandIndex(lab.diabetes, pugmm.diabetes.Gopt$label)}
# GPCMs
if(mclust.diabetes$G != length(unique(lab.diabetes))){
  mclust.diabetes.Gopt <- Mclust(x, G = length(unique(lab.diabetes)), control = emControl(itmax = 500))
  mclust.diabetes.Gopt$modelName
  nMclustParams(mclust.diabetes.Gopt$modelName, mclust.diabetes.Gopt$d, mclust.diabetes.Gopt$G)
  mclust.diabetes.Gopt.ari <- adjustedRandIndex(lab.diabetes, mclust.diabetes.Gopt$classification)}
# PGMMs
if(pgmm.diabetes$g != length(unique(lab.diabetes))){
  pgmm.diabetes.Gopt <- pgmmEM(x, rG = length(unique(lab.diabetes)), rq = 1:m.max, relax = TRUE)
  pgmm.diabetes.Gopt$q
  pgmm.diabetes.Gopt$model
  PGMM_dfree(pgmm.diabetes.Gopt$q, mclust.diabetes$d, pgmm.diabetes.Gopt$g, pgmm.diabetes.Gopt$model) 
  pgmm.diabetes.Gopt.ari <- adjustedRandIndex(lab.diabetes, pgmm.diabetes.Gopt$map)}
# HDDC
if(hddc.diabetes$K != length(unique(lab.diabetes))){
  hddc.diabetes.Gopt <- hddc(x, K = length(unique(lab.diabetes)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.diabetes.Gopt$model
  hddc.diabetes.Gopt$complexity
  hddc.diabetes.Gopt.ari <- adjustedRandIndex(lab.diabetes, hddc.diabetes.Gopt$class)}

################################ AIS ###########################################
data(ais, package = "sn")
x <- scale(ais[, -c(1,2)])
lab.ais <- as.numeric(ais[, 1])
G.max <- length(unique(lab.ais)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.ais <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.ais <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.ais <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.ais <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.ais$G
pugmm.ais$m
pugmm.ais$model.name
pugmm.ais$pm.free
pugmm.ais.ari <- adjustedRandIndex(lab.ais, pugmm.ais$label)
# GPCMs
mclust.ais$G
mclust.ais$modelName
nMclustParams(mclust.ais$modelName, mclust.ais$d, mclust.ais$G)
mclust.ais.ari <- adjustedRandIndex(lab.ais, mclust.ais$classification)
# PGMMs
pgmm.ais$g
pgmm.ais$q
pgmm.ais$model
PGMM_dfree(pgmm.ais$q, mclust.ais$d, pgmm.ais$g, pgmm.ais$model)
pgmm.ais.ari <- adjustedRandIndex(lab.ais, pgmm.ais$map)
# HDDC
hddc.ais$K
hddc.ais$d
hddc.ais$model
hddc.ais$complexity
hddc.ais.ari <- adjustedRandIndex(lab.ais, hddc.ais$class)

# Gopt
# PUGMMs
if(pugmm.ais$G != length(unique(lab.ais))){
  pugmm.ais.Gopt <- pugmm(x, G = length(unique(lab.ais)), m = 1:m.max)
  pugmm.ais.Gopt$m
  pugmm.ais.Gopt$model.name
  pugmm.ais.Gopt$pm.free
  pugmm.ais.Gopt.ari <- adjustedRandIndex(lab.ais, pugmm.ais.Gopt$label)}
# GPCMs
if(mclust.ais$G != length(unique(lab.ais))){
  mclust.ais.Gopt <- Mclust(x, G = length(unique(lab.ais)), control = emControl(itmax = 500))
  mclust.ais.Gopt$modelName
  nMclustParams(mclust.ais.Gopt$modelName, mclust.ais.Gopt$d, mclust.ais.Gopt$G)
  mclust.ais.Gopt.ari <- adjustedRandIndex(lab.ais, mclust.ais.Gopt$classification)}
# PGMMs
if(pgmm.ais$g != length(unique(lab.ais))){
  pgmm.ais.Gopt <- pgmmEM(x, rG = length(unique(lab.ais)), rq = 1:m.max, relax = TRUE)
  pgmm.ais.Gopt$q
  pgmm.ais.Gopt$model
  PGMM_dfree(pgmm.ais.Gopt$q, mclust.ais$d, pgmm.ais.Gopt$g, pgmm.ais.Gopt$model) 
  pgmm.ais.Gopt.ari <- adjustedRandIndex(lab.ais, pgmm.ais.Gopt$map)}
# HDDC
if(hddc.ais$K != length(unique(lab.ais))){
  hddc.ais.Gopt <- hddc(x, K = length(unique(lab.ais)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.ais.Gopt$model
  hddc.ais.Gopt$complexity
  hddc.ais.Gopt.ari <- adjustedRandIndex(lab.ais, hddc.ais.Gopt$class)}

################################ CERAMIC #######################################
ceramic <- read.csv("DIRECTORY_NAME/Ceramic.txt", header = FALSE)
x <- scale(ceramic[, -c(1,2)])
lab.ceramic <- ceramic[, 2]
G.max <- length(unique(lab.ceramic)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.ceramic <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.ceramic <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.ceramic <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.ceramic <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.ceramic$G
pugmm.ceramic$m
pugmm.ceramic$model.name
pugmm.ceramic$pm.free
pugmm.ceramic.ari <- adjustedRandIndex(lab.ceramic, pugmm.ceramic$label)
# GPCMs
mclust.ceramic$G
mclust.ceramic$modelName
nMclustParams(mclust.ceramic$modelName, mclust.ceramic$d, mclust.ceramic$G)
mclust.ceramic.ari <- adjustedRandIndex(lab.ceramic, mclust.ceramic$classification)
# PGMMs
pgmm.ceramic$g
pgmm.ceramic$q
pgmm.ceramic$model
PGMM_dfree(pgmm.ceramic$q, mclust.ceramic$d, pgmm.ceramic$g, pgmm.ceramic$model)
pgmm.ceramic.ari <- adjustedRandIndex(lab.ceramic, pgmm.ceramic$map)
# HDDC
hddc.ceramic$K
hddc.ceramic$d
hddc.ceramic$model
hddc.ceramic$complexity
hddc.ceramic.ari <- adjustedRandIndex(lab.ceramic, hddc.ceramic$class)

# Gopt
# PUGMMs
if(pugmm.ceramic$G != length(unique(lab.ceramic))){
  pugmm.ceramic.Gopt <- pugmm(x, G = length(unique(lab.ceramic)), m = 1:m.max)
  pugmm.ceramic.Gopt$m
  pugmm.ceramic.Gopt$model.name
  pugmm.ceramic.Gopt$pm.free
  pugmm.ceramic.Gopt.ari <- adjustedRandIndex(lab.ceramic, pugmm.ceramic.Gopt$label)}
# GPCMs
if(mclust.ceramic$G != length(unique(lab.ceramic))){
  mclust.ceramic.Gopt <- Mclust(x, G = length(unique(lab.ceramic)), control = emControl(itmax = 500))
  mclust.ceramic.Gopt$modelName
  nMclustParams(mclust.ceramic.Gopt$modelName, mclust.ceramic.Gopt$d, mclust.ceramic.Gopt$G)
  mclust.ceramic.Gopt.ari <- adjustedRandIndex(lab.ceramic, mclust.ceramic.Gopt$classification)}
# PGMMs
if(pgmm.ceramic$g != length(unique(lab.ceramic))){
  pgmm.ceramic.Gopt <- pgmmEM(x, rG = length(unique(lab.ceramic)), rq = 1:m.max, relax = TRUE)
  pgmm.ceramic.Gopt$q
  pgmm.ceramic.Gopt$model
  PGMM_dfree(pgmm.ceramic.Gopt$q, mclust.ceramic$d, pgmm.ceramic.Gopt$g, pgmm.ceramic.Gopt$model) 
  pgmm.ceramic.Gopt.ari <- adjustedRandIndex(lab.ceramic, pgmm.ceramic.Gopt$map)}
# HDDC
if(hddc.ceramic$K != length(unique(lab.ceramic))){
  hddc.ceramic.Gopt <- hddc(x, K = length(unique(lab.ceramic)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.ceramic.Gopt$model
  hddc.ceramic.Gopt$complexity
  hddc.ceramic.Gopt.ari <- adjustedRandIndex(lab.ceramic, hddc.ceramic.Gopt$class)}

################################ BANKNOTES #####################################
data(banknote, package = "mclust")
x <- scale(banknote[, -1])
lab.banknote <- as.numeric(banknote[, 1])
G.max <- length(unique(lab.banknote)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.banknote <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.banknote <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.banknote <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.banknote <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.banknote$G
pugmm.banknote$m
pugmm.banknote$model.name
pugmm.banknote$pm.free
pugmm.banknote.ari <- adjustedRandIndex(lab.banknote, pugmm.banknote$label)
# GPCMs
mclust.banknote$G
mclust.banknote$modelName
nMclustParams(mclust.banknote$modelName, mclust.banknote$d, mclust.banknote$G)
mclust.banknote.ari <- adjustedRandIndex(lab.banknote, mclust.banknote$classification)
# PGMMs
pgmm.banknote$g
pgmm.banknote$q
pgmm.banknote$model
PGMM_dfree(pgmm.banknote$q, mclust.banknote$d, pgmm.banknote$g, pgmm.banknote$model)
pgmm.banknote.ari <- adjustedRandIndex(lab.banknote, pgmm.banknote$map)
# HDDC
hddc.banknote$K
hddc.banknote$d
hddc.banknote$model
hddc.banknote$complexity
hddc.banknote.ari <- adjustedRandIndex(lab.banknote, hddc.banknote$class)

# Gopt
# PUGMMs
if(pugmm.banknote$G != length(unique(lab.banknote))){
  pugmm.banknote.Gopt <- pugmm(x, G = length(unique(lab.banknote)), m = 1:m.max)
  pugmm.banknote.Gopt$m
  pugmm.banknote.Gopt$model.name
  pugmm.banknote.Gopt$pm.free
  pugmm.banknote.Gopt.ari <- adjustedRandIndex(lab.banknote, pugmm.banknote.Gopt$label)}
# GPCMs
if(mclust.banknote$G != length(unique(lab.banknote))){
  mclust.banknote.Gopt <- Mclust(x, G = length(unique(lab.banknote)), control = emControl(itmax = 500))
  mclust.banknote.Gopt$modelName
  nMclustParams(mclust.banknote.Gopt$modelName, mclust.banknote.Gopt$d, mclust.banknote.Gopt$G)
  mclust.banknote.Gopt.ari <- adjustedRandIndex(lab.banknote, mclust.banknote.Gopt$classification)}
# PGMMs
if(pgmm.banknote$g != length(unique(lab.banknote))){
  pgmm.banknote.Gopt <- pgmmEM(x, rG = length(unique(lab.banknote)), rq = 1:m.max, relax = TRUE)
  pgmm.banknote.Gopt$q
  pgmm.banknote.Gopt$model
  PGMM_dfree(pgmm.banknote.Gopt$q, mclust.banknote$d, pgmm.banknote.Gopt$g, pgmm.banknote.Gopt$model) 
  pgmm.banknote.Gopt.ari <- adjustedRandIndex(lab.banknote, pgmm.banknote.Gopt$map)}
# HDDC
if(hddc.banknote$K != length(unique(lab.banknote))){
  hddc.banknote.Gopt <- hddc(x, K = length(unique(lab.banknote)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.banknote.Gopt$model
  hddc.banknote.Gopt$complexity
  hddc.banknote.Gopt.ari <- adjustedRandIndex(lab.banknote, hddc.banknote.Gopt$class)}

################################ COFFEE ########################################
data(coffee, package = "pgmm")
x <- scale(coffee[, -c(1, 2)])
lab.coffee <- coffee[, 1]
G.max <- length(unique(lab.coffee)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.coffee <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.coffee <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.coffee <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.coffee <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.coffee$G
pugmm.coffee$m
pugmm.coffee$model.name
pugmm.coffee$pm.free
pugmm.coffee.ari <- adjustedRandIndex(lab.coffee, pugmm.coffee$label)
# GPCMs
mclust.coffee$G
mclust.coffee$modelName
nMclustParams(mclust.coffee$modelName, mclust.coffee$d, mclust.coffee$G)
mclust.coffee.ari <- adjustedRandIndex(lab.coffee, mclust.coffee$classification)
# PGMMs
# pgmm.coffee$g
# pgmm.coffee$q
# pgmm.coffee$model
# PGMM_dfree(pgmm.coffee$q, mclust.coffee$d, pgmm.coffee$g, pgmm.coffee$model)
# pgmm.coffee.ari <- adjustedRandIndex(lab.coffee, pgmm.coffee$map)
# HDDC
hddc.coffee$K
hddc.coffee$d
hddc.coffee$model
hddc.coffee$complexity
hddc.coffee.ari <- adjustedRandIndex(lab.coffee, hddc.coffee$class)

# Gopt
# PUGMMs
if(pugmm.coffee$G != length(unique(lab.coffee))){
  pugmm.coffee.Gopt <- pugmm(x, G = length(unique(lab.coffee)), m = 1:m.max)
  pugmm.coffee.Gopt$m
  pugmm.coffee.Gopt$model.name
  pugmm.coffee.Gopt$pm.free
  pugmm.coffee.Gopt.ari <- adjustedRandIndex(lab.coffee, pugmm.coffee.Gopt$label)}
# GPCMs
if(mclust.coffee$G != length(unique(lab.coffee))){
  mclust.coffee.Gopt <- Mclust(x, G = length(unique(lab.coffee)), control = emControl(itmax = 500))
  mclust.coffee.Gopt$modelName
  nMclustParams(mclust.coffee.Gopt$modelName, mclust.coffee.Gopt$d, mclust.coffee.Gopt$G)
  mclust.coffee.Gopt.ari <- adjustedRandIndex(lab.coffee, mclust.coffee.Gopt$classification)}
# PGMMs
pgmm.coffee.Gopt <- pgmmEM(x, rG = length(unique(lab.coffee)), rq = 1:m.max, relax = TRUE)
pgmm.coffee.Gopt$q
pgmm.coffee.Gopt$model
PGMM_dfree(pgmm.coffee.Gopt$q, mclust.coffee$d, pgmm.coffee.Gopt$g, pgmm.coffee.Gopt$model) 
pgmm.coffee.Gopt.ari <- adjustedRandIndex(lab.coffee, pgmm.coffee.Gopt$map)
# HDDC
if(hddc.coffee$K != length(unique(lab.coffee))){
  hddc.coffee.Gopt <- hddc(x, K = length(unique(lab.coffee)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.coffee.Gopt$model
  hddc.coffee.Gopt$complexity
  hddc.coffee.Gopt.ari <- adjustedRandIndex(lab.coffee, hddc.coffee.Gopt$class)}

################################ SOBAR #########################################
sobar <- read.csv("sobar.csv", sep = ",")
x <- scale(sobar[, -20])
lab.sobar <- sobar[, 20]
G.max <- length(unique(lab.sobar)) + 2
if (dim(x)[2] <= 10){m.max <- dim(x)[2]} else {m.max <- 10}

pugmm.sobar <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.sobar <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.sobar <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.sobar <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.sobar$G
pugmm.sobar$m
pugmm.sobar$model.name
pugmm.sobar$pm.free
pugmm.sobar.ari <- adjustedRandIndex(lab.sobar, pugmm.sobar$label)
# GPCMs
mclust.sobar$G
mclust.sobar$modelName
nMclustParams(mclust.sobar$modelName, mclust.sobar$d, mclust.sobar$G)
mclust.sobar.ari <- adjustedRandIndex(lab.sobar, mclust.sobar$classification)
# PGMMs
pgmm.sobar$g
pgmm.sobar$q
pgmm.sobar$model
PGMM_dfree(pgmm.sobar$q, mclust.sobar$d, pgmm.sobar$g, pgmm.sobar$model)
pgmm.sobar.ari <- adjustedRandIndex(lab.sobar, pgmm.sobar$map)
# HDDC
hddc.sobar$K
hddc.sobar$d
hddc.sobar$model
hddc.sobar$complexity
hddc.sobar.ari <- adjustedRandIndex(lab.sobar, hddc.sobar$class)

# Gopt
# PUGMMs
if(pugmm.sobar$G != length(unique(lab.sobar))){
  pugmm.sobar.Gopt <- pugmm(x, G = length(unique(lab.sobar)), m = 1:m.max)
  pugmm.sobar.Gopt$m
  pugmm.sobar.Gopt$model.name
  pugmm.sobar.Gopt$pm.free
  pugmm.sobar.Gopt.ari <- adjustedRandIndex(lab.sobar, pugmm.sobar.Gopt$label)}
# GPCMs
if(mclust.sobar$G != length(unique(lab.sobar))){
  mclust.sobar.Gopt <- Mclust(x, G = length(unique(lab.sobar)), control = emControl(itmax = 500))
  mclust.sobar.Gopt$modelName
  nMclustParams(mclust.sobar.Gopt$modelName, mclust.sobar.Gopt$d, mclust.sobar.Gopt$G)
  mclust.sobar.Gopt.ari <- adjustedRandIndex(lab.sobar, mclust.sobar.Gopt$classification)}
# PGMMs
if(pgmm.sobar$g != length(unique(lab.sobar))){
  pgmm.sobar.Gopt <- pgmmEM(x, rG = length(unique(lab.sobar)), rq = 1:m.max, relax = TRUE)
  pgmm.sobar.Gopt$q
  pgmm.sobar.Gopt$model
  PGMM_dfree(pgmm.sobar.Gopt$q, mclust.sobar$d, pgmm.sobar.Gopt$g, pgmm.sobar.Gopt$model) 
  pgmm.sobar.Gopt.ari <- adjustedRandIndex(lab.sobar, pgmm.sobar.Gopt$map)}
# HDDC
if(hddc.sobar$K != length(unique(lab.sobar))){
  hddc.sobar.Gopt <- hddc(x, K = length(unique(lab.sobar)), model = "ALL", itermax = 500, mc.cores = 10, kmeans.control = list(nstart = 100))
  hddc.sobar.Gopt$model
  hddc.sobar.Gopt$complexity
  hddc.sobar.Gopt.ari <- adjustedRandIndex(lab.sobar, hddc.sobar.Gopt$class)}

