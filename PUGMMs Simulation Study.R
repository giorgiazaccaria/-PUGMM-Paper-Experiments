################################################################################
#                      Script for simulation study                             #
#                                                                              #
# In: Cavicchia, C. Vichi, M., Zaccaria, G. (2024+)                            #
# Parsimonious Gaussian mixture models. Revised version.                       #                                
################################################################################
# If necessary
# install.packages("devtools")
devtools:: install_github("PUGMM-authors/PUGMM")
library(PUGMM)
# If necessary
# install.packages("mclust")
library(mclust)
# If necessary
# install.packages("MASS")
library(MASS)
# If necessary
# install.packages("pracma")
library(pracma)
# If necessary
# install.packages("doParallel")
library(doParallel)
# If necessary
# install.packages("foreach")
library(foreach)
# Use the package "MixSim" for computing the overlapping level among components.
# If necessary
# install.packages("MixSim")
# library(MixSim)

################################################################################
##################
#    1.EUUU      #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
sigma.v <- 1.5
sigma.w <- 1.3
sigma.b <- 0.6
sigma.EUUU <- vector(mode = "list", length = param[3])
sigma.EUUU[1:param[3]] <- list(V.th[[1]] %*% (sigma.w * diag(param[4]) + sigma.b * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% ((sigma.v - sigma.w) * diag(param[4])) %*% t(V.th[[1]]))))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EUUU[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EUUU[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EUUU[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EUUU.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUUU.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUUU.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUUU.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EUUU.well[[s]]$res$label)
  if (pugmm.EUUU.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUUU.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUUU.well[[s]]$res$model.name == "EUUU") {
    count.model <- count.model + 1
  }
  stat.EUUU.well[s, 2] <- pugmm.EUUU.well[[s]]$runtime[3]
}
mstat.EUUU.well <- data.frame(count.G, mean(stat.EUUU.well[, 1]), count.m, count.model, mean(stat.EUUU.well[, 2]))
colnames(mstat.EUUU.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EUUU[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EUUU[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EUUU[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EUUU.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUUU.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUUU.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUUU.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EUUU.close[[s]]$res$label)
  if (pugmm.EUUU.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUUU.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUUU.close[[s]]$res$model.name == "EUUU") {
    count.model <- count.model + 1
  }
  stat.EUUU.close[s, 2] <- pugmm.EUUU.close[[s]]$runtime[3]
}
mstat.EUUU.close <- data.frame(count.G, mean(stat.EUUU.close[, 1]), count.m, count.model, mean(stat.EUUU.close[, 2]))
colnames(mstat.EUUU.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     2.EUUE     #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
sigma.v <- 1.5
sigma.w <- 1.3
Sigma.b <- matrix(0, param[4], param[4])
Sigma.b[2, 3] = Sigma.b[3, 2] <- 0.9
Sigma.b[1, 2] = Sigma.b[1, 3] = Sigma.b[2, 1] = Sigma.b[3, 1] <- 0.5
sigma.EUUE <- vector(mode = "list", length = param[3])
sigma.EUUE[1:param[3]] <- list(V.th[[1]] %*% (sigma.w * diag(param[4]) + Sigma.b) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% ((sigma.v - sigma.w) * diag(param[4])) %*% t(V.th[[1]]))))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EUUE[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EUUE[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EUUE[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EUUE.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUUE.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUUE.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUUE.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EUUE.well[[s]]$res$label)
  if (pugmm.EUUE.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUUE.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUUE.well[[s]]$res$model.name == "EUUE") {
    count.model <- count.model + 1
  }
  stat.EUUE.well[s, 2] <- pugmm.EUUE.well[[s]]$runtime[3]
}
mstat.EUUE.well <- data.frame(count.G, mean(stat.EUUE.well[, 1]), count.m, count.model, mean(stat.EUUE.well[, 2]))
colnames(mstat.EUUE.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EUUE[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EUUE[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EUUE[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EUUE.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUUE.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUUE.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUUE.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EUUE.close[[s]]$res$label)
  if (pugmm.EUUE.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUUE.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUUE.close[[s]]$res$model.name == "EUUE") {
    count.model <- count.model + 1
  }
  stat.EUUE.close[s, 2] <- pugmm.EUUE.close[[s]]$runtime[3]
}
mstat.EUUE.close <- data.frame(count.G, mean(stat.EUUE.close[, 1]), count.m, count.model, mean(stat.EUUE.close[, 2]))
colnames(mstat.EUUE.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     3.EUEE     #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
sigma.v <- 2
Sigma.w <- matrix(0,  param[4], param[4])
Sigma.w[1, 1] <- 1.2
Sigma.w[2, 2] <- 0.5
Sigma.w[3, 3] <- 1.4
Sigma.b <- matrix(0, param[4], param[4])
Sigma.b <- matrix(0, param[4], param[4])
Sigma.b[2, 3] = Sigma.b[3, 2] <- 0.5
Sigma.b[1, 2] = Sigma.b[1, 3] = Sigma.b[2, 1] = Sigma.b[3, 1] <- 0
sigma.EUEE <- vector(mode = "list", length = param[3])
sigma.EUEE[1:param[3]] <- list(V.th[[1]] %*% (Sigma.w + Sigma.b) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (sigma.v * diag(param[4]) - Sigma.w)  %*% t(V.th[[1]]))))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EUEE[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EUEE[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EUEE[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EUEE.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUEE.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUEE.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUEE.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EUEE.well[[s]]$res$label)
  if (pugmm.EUEE.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUEE.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUEE.well[[s]]$res$model.name == "EUEE") {
    count.model <- count.model + 1
  }
  stat.EUEE.well[s, 2] <- pugmm.EUEE.well[[s]]$runtime[3]
}
mstat.EUEE.well <- data.frame(count.G, mean(stat.EUEE.well[, 1]), count.m, count.model, mean(stat.EUEE.well[, 2]))
colnames(mstat.EUEE.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(2.8, param[2]), rep(6.8, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EUEE[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EUEE[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EUEE[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EUEE.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EUEE.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EUEE.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EUEE.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EUEE.close[[s]]$res$label)
  if (pugmm.EUEE.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EUEE.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EUEE.close[[s]]$res$model.name == "EUEE") {
    count.model <- count.model + 1
  }
  stat.EUEE.close[s, 2] <- pugmm.EUEE.close[[s]]$runtime[3]
}
mstat.EUEE.close <- data.frame(count.G, mean(stat.EUEE.close[, 1]), count.m, count.model, mean(stat.EUEE.close[, 2]))
colnames(mstat.EUEE.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     4.EEEU     #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
Sigma.v <- matrix(0, param[4], param[4])
Sigma.v[1, 1] <- 1.7
Sigma.v[2, 2] <- 1.5
Sigma.v[3, 3] <- 1.4
Sigma.w <- matrix(0,  param[4], param[4])
Sigma.w[1, 1] <- 1.5
Sigma.w[2, 2] <- 1.1
Sigma.w[3, 3] <- 1.2
sigma.b <- 0.60
sigma.EEEU <- vector(mode = "list", length = param[3])
sigma.EEEU[1:param[3]] <- list(V.th[[1]] %*% (Sigma.w + sigma.b * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v - Sigma.w)  %*% t(V.th[[1]]))))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EEEU[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EEEU[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EEEU[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EEEU.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEU.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(X.well[[i]], G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEU.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEU.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEU.well[[s]]$res$label)
  if (pugmm.EEEU.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEU.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEU.well[[s]]$res$model.name == "EEEU") {
    count.model <- count.model + 1
  }
  stat.EEEU.well[s, 2] <- pugmm.EEEU.well[[s]]$runtime[3]
}
mstat.EEEU.well <- data.frame(count.G, mean(stat.EEEU.well[, 1]), count.m, count.model, mean(stat.EEEU.well[, 2]))
colnames(mstat.EEEU.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EEEU[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EEEU[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EEEU[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EEEU.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEU.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEU.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEU.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEU.close[[s]]$res$label)
  if (pugmm.EEEU.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEU.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEU.close[[s]]$res$model.name == "EEEU") {
    count.model <- count.model + 1
  }
  stat.EEEU.close[s, 2] <- pugmm.EEEU.close[[s]]$runtime[3]
}
mstat.EEEU.close <- data.frame(count.G, mean(stat.EEEU.close[, 1]), count.m, count.model, mean(stat.EEEU.close[, 2]))
colnames(mstat.EEEU.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     5.EEEE     #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
Sigma.v <- matrix(0, param[4], param[4])
Sigma.v[1, 1] <- 1.8
Sigma.v[2, 2] <- 1.7
Sigma.v[3, 3] <- 1.6
Sigma.w <- matrix(0,  param[4], param[4])
Sigma.w[1, 1] <- 1.7
Sigma.w[2, 2] <- 1.3
Sigma.w[3, 3] <- 1.4
Sigma.b <- matrix(0, param[4], param[4])
Sigma.b <- matrix(0, param[4], param[4])
Sigma.b[2, 3] = Sigma.b[3, 2] <- 0.9
Sigma.b[1, 2] = Sigma.b[1, 3] = Sigma.b[2, 1] = Sigma.b[3, 1] <- 0.3
sigma.EEEE <- vector(mode = "list", length = param[3])
sigma.EEEE[1:param[3]] <- list(V.th[[1]] %*% (Sigma.w + Sigma.b) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v - Sigma.w)  %*% t(V.th[[1]]))))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EEEE[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EEEE[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EEEE[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EEEE.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEE.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEE.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEE.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEE.well[[s]]$res$label)
  if (pugmm.EEEE.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEE.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEE.well[[s]]$res$model.name == "EEEE") {
    count.model <- count.model + 1
  }
  stat.EEEE.well[s, 2] <- pugmm.EEEE.well[[s]]$runtime[3]
}
mstat.EEEE.well <- data.frame(count.G, mean(stat.EEEE.well[, 1]), count.m, count.model, mean(stat.EEEE.well[, 2]))
colnames(mstat.EEEE.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EEEE[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EEEE[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EEEE[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EEEE.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEE.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEE.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEE.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEE.close[[s]]$res$label)
  if (pugmm.EEEE.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEE.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEE.close[[s]]$res$model.name == "EEEE") {
    count.model <- count.model + 1
  }
  stat.EEEE.close[s, 2] <- pugmm.EEEE.close[[s]]$runtime[3]
}
mstat.EEEE.close <- data.frame(count.G, mean(stat.EEEE.close[, 1]), count.m, count.model, mean(stat.EEEE.close[, 2]))
colnames(mstat.EEEE.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     6.EEEF    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
Sigma.v <- vector(mode = "list", length = param[3])
Sigma.v <- matrix(0, param[4], param[4])
Sigma.v[1, 1] <- 1.8
Sigma.v[2, 2] <- 1.7
Sigma.v[3, 3] <- 1.6
Sigma.w <- vector(mode = "list", length = param[3])
Sigma.w[[1]] <- matrix(0,  param[4], param[4])
Sigma.w[[1]][1, 1] <- 1.7
Sigma.w[[1]][2, 2] <- 1.3
Sigma.w[[1]][3, 3] <- 1.4
Sigma.w[[2]] <- matrix(0,  param[4], param[4])
Sigma.w[[2]][1, 1] <- 1.7
Sigma.w[[2]][2, 2] <- 1.3
Sigma.w[[2]][3, 3] <- 1.4
Sigma.w[[3]] <- matrix(0,  param[4], param[4])
Sigma.w[[3]][1, 1] <- 1.7
Sigma.w[[3]][2, 2] <- 1.3
Sigma.w[[3]][3, 3] <- 1.4
Sigma.b <- vector(mode = "list", length = param[3])
Sigma.b[[1]] <- matrix(0, param[4], param[4])
Sigma.b[[1]][2, 3] = Sigma.b[[1]][3, 2] <- 0.8
Sigma.b[[1]][1, 2] = Sigma.b[[1]][1, 3] = Sigma.b[[1]][2, 1] = Sigma.b[[1]][3, 1] <- 0.3
Sigma.b[[2]] <- matrix(0, param[4], param[4])
Sigma.b[[2]][2, 3] = Sigma.b[[2]][3, 2] <- 1.3
Sigma.b[[2]][1, 2] = Sigma.b[[2]][1, 3] = Sigma.b[[2]][2, 1] = Sigma.b[[2]][3, 1] <- -0.2
Sigma.b[[3]] <- matrix(0, param[4], param[4])
Sigma.b[[3]][2, 3] = Sigma.b[[3]][3, 2] <- 0.7
Sigma.b[[3]][1, 2] = Sigma.b[[3]][1, 3] = Sigma.b[[3]][2, 1] = Sigma.b[[3]][3, 1] <- 0
# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EEEF[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EEEF[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EEEF[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EEEF.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEF.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEF.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEF.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEF.well[[s]]$res$label)
  if (pugmm.EEEF.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEF.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEF.well[[s]]$res$model.name == "EEEF") {
    count.model <- count.model + 1
  }
  stat.EEEF.well[s, 2] <- pugmm.EEEF.well[[s]]$runtime[3]
}
mstat.EEEF.well <- data.frame(count.G, mean(stat.EEEF.well[, 1]), count.m, count.model, mean(stat.EEEF.well[, 2]))
colnames(mstat.EEEF.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EEEF[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EEEF[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EEEF[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EEEF.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEEF.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEEF.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEEF.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EEEF.close[[s]]$res$label)
  if (pugmm.EEEF.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEEF.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEEF.close[[s]]$res$model.name == "EEEF") {
    count.model <- count.model + 1
  }
  stat.EEEF.close[s, 2] <- pugmm.EEEF.close[[s]]$runtime[3]
}
mstat.EEEF.close <- data.frame(count.G, mean(stat.EEEF.close[, 1]), count.m, count.model, mean(stat.EEEF.close[, 2]))
colnames(mstat.EEEF.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     7.EEFF    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
Sigma.v <- matrix(0, param[4], param[4])
Sigma.v[1, 1] <- 2.9
Sigma.v[2, 2] <- 3.0
Sigma.v[3, 3] <- 3.1
Sigma.w <- vector(mode = "list", length = param[3])
Sigma.w[[1]] <- matrix(0,  param[4], param[4])
Sigma.w[[1]][1, 1] <- 1.9
Sigma.w[[1]][2, 2] <- 0.4
Sigma.w[[1]][3, 3] <- 1.4
Sigma.w[[2]] <- matrix(0,  param[4], param[4])
Sigma.w[[2]][1, 1] <- 0.2
Sigma.w[[2]][2, 2] <- 1.5
Sigma.w[[2]][3, 3] <- 2.4
Sigma.w[[3]] <- matrix(0,  param[4], param[4])
Sigma.w[[3]][1, 1] <- 0.2
Sigma.w[[3]][2, 2] <- 1.8
Sigma.w[[3]][3, 3] <- 0.8
Sigma.b <- vector(mode = "list", length = param[3])
Sigma.b[[1]] <- matrix(0, param[4], param[4])
Sigma.b[[1]][2, 3] = Sigma.b[[1]][3, 2] <- 0.3
Sigma.b[[1]][1, 2] = Sigma.b[[1]][1, 3] = Sigma.b[[1]][2, 1] = Sigma.b[[1]][3, 1] <- -0.5
Sigma.b[[2]] <- matrix(0, param[4], param[4])
Sigma.b[[2]][2, 3] = Sigma.b[[2]][3, 2] <- 0
Sigma.b[[2]][1, 2] = Sigma.b[[2]][1, 3] = Sigma.b[[2]][2, 1] = Sigma.b[[2]][3, 1] <- -0.3
Sigma.b[[3]] <- matrix(0, param[4], param[4])
Sigma.b[[3]][2, 3] = Sigma.b[[3]][3, 2] <- 0.2
Sigma.b[[3]][1, 2] = Sigma.b[[3]][1, 3] = Sigma.b[[3]][2, 1] = Sigma.b[[3]][3, 1] <- 0
sigma.EEFF <- vector(mode = "list", length = param[3])
sigma.EEFF[[1]] <- V.th[[1]] %*% (Sigma.w[[1]] + Sigma.b[[1]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v - Sigma.w[[1]])  %*% t(V.th[[1]])))
sigma.EEFF[[2]] <- V.th[[1]] %*% (Sigma.w[[2]] + Sigma.b[[2]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v - Sigma.w[[2]])  %*% t(V.th[[1]])))
sigma.EEFF[[3]] <- V.th[[1]] %*% (Sigma.w[[3]] + Sigma.b[[3]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v - Sigma.w[[3]])  %*% t(V.th[[1]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EEFF[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EEFF[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EEFF[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EEFF.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEFF.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEFF.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEFF.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EEFF.well[[s]]$res$label)
  if (pugmm.EEFF.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEFF.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEFF.well[[s]]$res$model.name == "EEFF") {
    count.model <- count.model + 1
  }
  stat.EEFF.well[s, 2] <- pugmm.EEFF.well[[s]]$runtime[3]
}
mstat.EEFF.well <- data.frame(count.G, mean(stat.EEFF.well[, 1]), count.m, count.model, mean(stat.EEFF.well[, 2]))
colnames(mstat.EEFF.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EEFF[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EEFF[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EEFF[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EEFF.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EEFF.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EEFF.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EEFF.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EEFF.close[[s]]$res$label)
  if (pugmm.EEFF.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EEFF.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EEFF.close[[s]]$res$model.name == "EEFF") {
    count.model <- count.model + 1
  }
  stat.EEFF.close[s, 2] <- pugmm.EEFF.close[[s]]$runtime[3]
}
mstat.EEFF.close <- data.frame(count.G, mean(stat.EEFF.close[, 1]), count.m, count.model, mean(stat.EEFF.close[, 2]))
colnames(mstat.EEFF.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     8.EFFF    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[1:param[3]] <- list(as.matrix(cbind(c(rep(1,3), rep(0,7)),
                                         c(rep(0,3), rep(1,3), rep(0,4)),
                                         c(rep(0,6), rep(1,4)))))
Sigma.v <- vector(mode = "list", length = param[3])
Sigma.v[[1]] <- matrix(0,  param[4], param[4])
Sigma.v[[1]][1, 1] <- 2.9
Sigma.v[[1]][2, 2] <- 3.0
Sigma.v[[1]][3, 3] <- 3.1
Sigma.v[[2]] <- matrix(0, param[4], param[4])
Sigma.v[[2]][1, 1] <- 2.7
Sigma.v[[2]][2, 2] <- 2.9
Sigma.v[[2]][3, 3] <- 2.5
Sigma.v[[3]] <- matrix(0, param[4], param[4])
Sigma.v[[3]][1, 1] <- 2.3
Sigma.v[[3]][2, 2] <- 2.0
Sigma.v[[3]][3, 3] <- 2.4
Sigma.w <- vector(mode = "list", length = param[3])
Sigma.w[[1]] <- matrix(0,  param[4], param[4])
Sigma.w[[1]][1, 1] <- 1.9
Sigma.w[[1]][2, 2] <- 0.4
Sigma.w[[1]][3, 3] <- 1.4
Sigma.w[[2]] <- matrix(0,  param[4], param[4])
Sigma.w[[2]][1, 1] <- 0.2
Sigma.w[[2]][2, 2] <- 1.5
Sigma.w[[2]][3, 3] <- 2.4
Sigma.w[[3]] <- matrix(0,  param[4], param[4])
Sigma.w[[3]][1, 1] <- 0.2
Sigma.w[[3]][2, 2] <- 1.8
Sigma.w[[3]][3, 3] <- 0.8
Sigma.b <- vector(mode = "list", length = param[3])
Sigma.b[[1]] <- matrix(0, param[4], param[4])
Sigma.b[[1]][2, 3] = Sigma.b[[1]][3, 2] <- 0.3
Sigma.b[[1]][1, 2] = Sigma.b[[1]][1, 3] = Sigma.b[[1]][2, 1] = Sigma.b[[1]][3, 1] <- -0.5
Sigma.b[[2]] <- matrix(0, param[4], param[4])
Sigma.b[[2]][2, 3] = Sigma.b[[2]][3, 2] <- 0
Sigma.b[[2]][1, 2] = Sigma.b[[2]][1, 3] = Sigma.b[[2]][2, 1] = Sigma.b[[2]][3, 1] <- -0.3
Sigma.b[[3]] <- matrix(0, param[4], param[4])
Sigma.b[[3]][2, 3] = Sigma.b[[3]][3, 2] <- 0.2
Sigma.b[[3]][1, 2] = Sigma.b[[3]][1, 3] = Sigma.b[[3]][2, 1] = Sigma.b[[3]][3, 1] <- 0
sigma.EFFF <- vector(mode = "list", length = param[3])
sigma.EFFF[[1]] <- V.th[[1]] %*% (Sigma.w[[1]] + Sigma.b[[1]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v[[1]] - Sigma.w[[1]])  %*% t(V.th[[1]])))
sigma.EFFF[[2]] <- V.th[[1]] %*% (Sigma.w[[2]] + Sigma.b[[2]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v[[2]] - Sigma.w[[2]])  %*% t(V.th[[1]])))
sigma.EFFF[[3]] <- V.th[[1]] %*% (Sigma.w[[3]] + Sigma.b[[3]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v[[3]] - Sigma.w[[3]])  %*% t(V.th[[1]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.EFFF[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.EFFF[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.EFFF[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.EFFF.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EFFF.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EFFF.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EFFF.well[s, 1] <- adjustedRandIndex(label.th, pugmm.EFFF.well[[s]]$res$label)
  if (pugmm.EFFF.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EFFF.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EFFF.well[[s]]$res$model.name == "EFFF") {
    count.model <- count.model + 1
  }
  stat.EFFF.well[s, 2] <- pugmm.EFFF.well[[s]]$runtime[3]
}
mstat.EFFF.well <- data.frame(count.G, mean(stat.EFFF.well[, 1]), count.m, count.model, mean(stat.EFFF.well[, 2]))
colnames(mstat.EFFF.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(1.8, param[2]), rep(3.2, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.EFFF[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.EFFF[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.EFFF[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.EFFF.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.EFFF.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.EFFF.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.EFFF.close[s, 1] <- adjustedRandIndex(label.th, pugmm.EFFF.close[[s]]$res$label)
  if (pugmm.EFFF.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.EFFF.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.EFFF.close[[s]]$res$model.name == "EFFF") {
    count.model <- count.model + 1
  }
  stat.EFFF.close[s, 2] <- pugmm.EFFF.close[[s]]$runtime[3]
}
mstat.EFFF.close <- data.frame(count.G, mean(stat.EFFF.close[, 1]), count.m, count.model, mean(stat.EFFF.close[, 2]))
colnames(mstat.EFFF.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#    9.FIII      #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[[1]] <- as.matrix(cbind(c(rep(1,3), rep(0,7)),
                             c(rep(0,3), rep(1,3), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[2]] <- as.matrix(cbind(c(rep(1,2), rep(0,8)),
                             c(rep(0,2), rep(1,4), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[3]] <- as.matrix(cbind(c(rep(1,4), rep(0,6)),
                             c(rep(0,4), rep(1,3), rep(0,3)),
                             c(rep(0,7), rep(1,3))))

sigma.v <- c(1.5, 2, 1.8)
sigma.w <- c(1.3, 1.7, 1.5)
sigma.b <- c(0.6, 0.9, 0.5)
sigma.FIII <- vector(mode = "list", length = param[3])
sigma.FIII[[1]] <- V.th[[1]] %*% (sigma.w[1] * diag(param[4]) + sigma.b[1] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% ((sigma.v[1] - sigma.w[1]) * diag(param[4])) %*% t(V.th[[1]])))
sigma.FIII[[2]] <- V.th[[2]] %*% (sigma.w[2] * diag(param[4]) + sigma.b[2] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[2]]) + diag(diag(V.th[[2]] %*% ((sigma.v[2] - sigma.w[2]) * diag(param[4])) %*% t(V.th[[2]])))
sigma.FIII[[3]] <- V.th[[3]] %*% (sigma.w[3] * diag(param[4]) + sigma.b[3] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[3]]) + diag(diag(V.th[[3]] %*% ((sigma.v[3] - sigma.w[3]) * diag(param[4])) %*% t(V.th[[3]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.FIII[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.FIII[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.FIII[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.FIII.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIII.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIII.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIII.well[s, 1] <- adjustedRandIndex(label.th, pugmm.FIII.well[[s]]$res$label)
  if (pugmm.FIII.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIII.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIII.well[[s]]$res$model.name == "FIII") {
    count.model <- count.model + 1
  }
  stat.FIII.well[s, 2] <- pugmm.FIII.well[[s]]$runtime[3]
}
mstat.FIII.well <- data.frame(count.G, mean(stat.FIII.well[, 1]), count.m, count.model, mean(stat.FIII.well[, 2]))
colnames(mstat.FIII.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3.5, param[2]), rep(7.5, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.FIII[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.FIII[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.FIII[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.FIII.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIII.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIII.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIII.close[s, 1] <- adjustedRandIndex(label.th, stat.FIII.close[[s]]$res$label)
  if (pugmm.FIII.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIII.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIII.close[[s]]$res$model.name == "FIII") {
    count.model <- count.model + 1
  }
  stat.FIII.close[s, 2] <- pugmm.FIII.close[[s]]$runtime[3]
}
mstat.FIII.close <- data.frame(count.G, mean(stat.FIII.close[, 1]), count.m, count.model, mean(stat.FIII.close[, 2]))
colnames(mstat.FIII.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     10.FIIF    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[[1]] <- as.matrix(cbind(c(rep(1,3), rep(0,7)),
                             c(rep(0,3), rep(1,3), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[2]] <- as.matrix(cbind(c(rep(1,2), rep(0,8)),
                             c(rep(0,2), rep(1,4), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[3]] <- as.matrix(cbind(c(rep(1,4), rep(0,6)),
                             c(rep(0,4), rep(1,3), rep(0,3)),
                             c(rep(0,7), rep(1,3))))
sigma.v <- c(3, 2.2, 4)
sigma.w <- c(2.6, 1.8, 3.5)
Sigma.b <- vector(mode = "list", length = param[3])
Sigma.b[[1]] <- matrix(0, param[4], param[4])
Sigma.b[[1]][2, 3] = Sigma.b[[1]][3, 2] <- 0.9
Sigma.b[[1]][1, 2] = Sigma.b[[1]][1, 3] = Sigma.b[[1]][2, 1] = Sigma.b[[1]][3, 1] <- 0.3
Sigma.b[[2]] <- matrix(0, param[4], param[4])
Sigma.b[[2]][2, 3] = Sigma.b[[2]][3, 2] <- 0.9
Sigma.b[[2]][1, 2] = Sigma.b[[2]][1, 3] = Sigma.b[[2]][2, 1] = Sigma.b[[2]][3, 1] <- -0.8
Sigma.b[[3]] <- matrix(0, param[4], param[4])
Sigma.b[[3]][2, 3] = Sigma.b[[3]][3, 2] <- 0.5
Sigma.b[[3]][1, 2] = Sigma.b[[3]][1, 3] = Sigma.b[[3]][2, 1] = Sigma.b[[3]][3, 1] <- 0
sigma.FIIF <- vector(mode = "list", length = param[3])
sigma.FIIF[[1]] <- V.th[[1]] %*% (sigma.w[1] * diag(param[4]) + Sigma.b[[1]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% ((sigma.v[1] - sigma.w[1]) * diag(param[4])) %*% t(V.th[[1]])))
sigma.FIIF[[2]] <- V.th[[2]] %*% (sigma.w[2] * diag(param[4]) + Sigma.b[[2]]) %*% t(V.th[[2]]) + diag(diag(V.th[[2]] %*% ((sigma.v[2] - sigma.w[2]) * diag(param[4])) %*% t(V.th[[2]])))
sigma.FIIF[[3]] <- V.th[[3]] %*% (sigma.w[3] * diag(param[4]) + Sigma.b[[3]]) %*% t(V.th[[3]]) + diag(diag(V.th[[3]] %*% ((sigma.v[3] - sigma.w[3]) * diag(param[4])) %*% t(V.th[[3]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.FIIF[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.FIIF[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.FIIF[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.FIIF.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIIF.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIIF.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIIF.well[s, 1] <- adjustedRandIndex(label.th, pugmm.FIIF.well[[s]]$res$label)
  if (pugmm.FIIF.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIIF.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIIF.well[[s]]$res$model.name == "FIIF") {
    count.model <- count.model + 1
  }
  stat.FIIF.well[s, 2] <- pugmm.FIIF.well[[s]]$runtime[3]
}
mstat.FIIF.well <- data.frame(count.G, mean(stat.FIIF.well[, 1]), count.m, count.model, mean(stat.FIIF.well[, 2]))
colnames(mstat.FIIF.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3, param[2]), rep(7, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.FIIF[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.FIIF[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.FIIF[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.FIIF.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIIF.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIIF.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIIF.close[s, 1] <- adjustedRandIndex(label.th, pugmm.FIIF.close[[s]]$res$label)
  if (pugmm.FIIF.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIIF.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIIF.close[[s]]$res$model.name == "FIIF") {
    count.model <- count.model + 1
  }
  stat.FIIF.close[s, 2] <- pugmm.FIIF.close[[s]]$runtime[3]
}
mstat.FIIF.close <- data.frame(count.G, mean(stat.FIIF.close[, 1]), count.m, count.model, mean(stat.FIIF.close[, 2]))
colnames(mstat.FIIF.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     11.FIFF    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[[1]] <- as.matrix(cbind(c(rep(1,3), rep(0,7)),
                             c(rep(0,3), rep(1,2), rep(0,5)),
                             c(rep(0,5), rep(1,5))))
V.th[[2]] <- as.matrix(cbind(c(rep(1,6), rep(0,4)),
                             c(rep(0,6), rep(1,2), rep(0,2)),
                             c(rep(0,8), rep(1,2))))
V.th[[3]] <- as.matrix(cbind(c(rep(1,2), rep(0,8)),
                             c(rep(0,2), rep(1,4), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
sigma.v <- c(2.7, 3.1, 3)
Sigma.w <- vector(mode = "list", length = param[3])
Sigma.w[[1]] <- matrix(0,  param[4], param[4])
Sigma.w[[1]][1, 1] <- 1.9
Sigma.w[[1]][2, 2] <- 0.4
Sigma.w[[1]][3, 3] <- 1.4
Sigma.w[[2]] <- matrix(0,  param[4], param[4])
Sigma.w[[2]][1, 1] <- 0.2
Sigma.w[[2]][2, 2] <- 1.5
Sigma.w[[2]][3, 3] <- 2.4
Sigma.w[[3]] <- matrix(0,  param[4], param[4])
Sigma.w[[3]][1, 1] <- 0.2
Sigma.w[[3]][2, 2] <- 1.8
Sigma.w[[3]][3, 3] <- 0.8
Sigma.b <- vector(mode = "list", length = param[3])
Sigma.b[[1]] <- matrix(0, param[4], param[4])
Sigma.b[[1]][2, 3] = Sigma.b[[1]][3, 2] <- 0.3
Sigma.b[[1]][1, 2] = Sigma.b[[1]][1, 3] = Sigma.b[[1]][2, 1] = Sigma.b[[1]][3, 1] <- -0.5
Sigma.b[[2]] <- matrix(0, param[4], param[4])
Sigma.b[[2]][2, 3] = Sigma.b[[2]][3, 2] <- 0
Sigma.b[[2]][1, 2] = Sigma.b[[2]][1, 3] = Sigma.b[[2]][2, 1] = Sigma.b[[2]][3, 1] <- -0.3
Sigma.b[[3]] <- matrix(0, param[4], param[4])
Sigma.b[[3]][2, 3] = Sigma.b[[3]][3, 2] <- 0.2
Sigma.b[[3]][1, 2] = Sigma.b[[3]][1, 3] = Sigma.b[[3]][2, 1] = Sigma.b[[3]][3, 1] <- 0
sigma.FIFF <- vector(mode = "list", length = param[3])
sigma.FIFF[[1]] <- V.th[[1]] %*% (Sigma.w[[1]] + Sigma.b[[1]]) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (sigma.v[1] * diag(param[4]) - Sigma.w[[1]]) %*% t(V.th[[1]])))
sigma.FIFF[[2]] <- V.th[[2]] %*% (Sigma.w[[2]] + Sigma.b[[2]]) %*% t(V.th[[2]]) + diag(diag(V.th[[2]] %*% (sigma.v[2] * diag(param[4]) - Sigma.w[[2]]) %*% t(V.th[[2]])))
sigma.FIFF[[3]] <- V.th[[3]] %*% (Sigma.w[[3]] + Sigma.b[[3]]) %*% t(V.th[[3]]) + diag(diag(V.th[[3]] %*% (sigma.v[3] * diag(param[4]) - Sigma.w[[3]]) %*% t(V.th[[3]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.FIFF[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.FIFF[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.FIFF[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.FIFF.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIFF.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIFF.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIFF.well[s, 1] <- adjustedRandIndex(label.th, pugmm.FIFF.well[[s]]$res$label)
  if (pugmm.FIFF.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIFF.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIFF.well[[s]]$res$model.name == "FIFF") {
    count.model <- count.model + 1
  }
  stat.FIFF.well[s, 2] <- pugmm.FIFF.well[[s]]$runtime[3]
}
mstat.FIFF.well <- data.frame(count.G, mean(stat.FIFF.well[, 1]), count.m, count.model, mean(stat.FIFF.well[, 2]))
colnames(mstat.FIFF.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(2.8, param[2]), rep(5.3, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.FIFF[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.FIFF[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.FIFF[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.FIFF.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FIFF.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FIFF.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FIFF.close[s, 1] <- adjustedRandIndex(label.th, pugmm.FIFF.close[[s]]$res$label)
  if (pugmm.FIFF.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FIFF.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FIFF.close[[s]]$res$model.name == "FIFF") {
    count.model <- count.model + 1
  }
  stat.FIFF.close[s, 2] <- pugmm.FIFF.close[[s]]$runtime[3]
}
mstat.FIFF.close <- data.frame(count.G, mean(stat.FIFF.close[, 1]), count.m, count.model, mean(stat.FIFF.close[, 2]))
colnames(mstat.FIFF.close) <- c("G", "mARI", "m", "model", "time")

################################################################################
##################
#     12.FFFI    #
##################
rm(list = ls())
cat("\014")

### Simulation study setting common to all cases
param <- c(n.obs = 200, p = 10, G = 3, m = 3, pp = c(0.25, 0.25, 0.5), nsamp = 100)
label.th <- c(rep(1, param[1]*param[5]), rep(2, param[1]*param[6]), rep(3, param[1]*param[7]))

### Covariance structure
V.th <- vector(mode = "list", length = param[3])
V.th[[1]] <- as.matrix(cbind(c(rep(1,3), rep(0,7)),
                             c(rep(0,3), rep(1,3), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[2]] <- as.matrix(cbind(c(rep(1,2), rep(0,8)),
                             c(rep(0,2), rep(1,4), rep(0,4)),
                             c(rep(0,6), rep(1,4))))
V.th[[3]] <- as.matrix(cbind(c(rep(1,4), rep(0,6)),
                             c(rep(0,4), rep(1,3), rep(0,3)),
                             c(rep(0,7), rep(1,3))))
Sigma.v <- vector(mode = "list", length = param[3])
Sigma.v[[1]] <- matrix(0,  param[4], param[4])
Sigma.v[[1]][1, 1] <- 3.3
Sigma.v[[1]][2, 2] <- 3
Sigma.v[[1]][3, 3] <- 2.7
Sigma.v[[2]] <- matrix(0,  param[4], param[4])
Sigma.v[[2]][1, 1] <- 2.5
Sigma.v[[2]][2, 2] <- 2.2
Sigma.v[[2]][3, 3] <- 1.9
Sigma.v[[3]] <- matrix(0,  param[4], param[4])
Sigma.v[[3]][1, 1] <- 4.3
Sigma.v[[3]][2, 2] <- 4
Sigma.v[[3]][3, 3] <- 3.7
Sigma.w <- vector(mode = "list", length = param[3])
Sigma.w[[1]] <- matrix(0,  param[4], param[4])
Sigma.w[[1]][1, 1] <- 2.9
Sigma.w[[1]][2, 2] <- 2.7
Sigma.w[[1]][3, 3] <- 2.4
Sigma.w[[2]] <- matrix(0,  param[4], param[4])
Sigma.w[[2]][1, 1] <- 1.8
Sigma.w[[2]][2, 2] <- 1.5
Sigma.w[[2]][3, 3] <- 1.2
Sigma.w[[3]] <- matrix(0,  param[4], param[4])
Sigma.w[[3]][1, 1] <- 3.8
Sigma.w[[3]][2, 2] <- 3.3
Sigma.w[[3]][3, 3] <- 3.5
sigma.b <- c(0.6, 0.9, 0.5)
sigma.FFFI <- vector(mode = "list", length = param[3])
sigma.FFFI[[1]] <- V.th[[1]] %*% (Sigma.w[[1]] + sigma.b[1] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[1]]) + diag(diag(V.th[[1]] %*% (Sigma.v[[1]] - Sigma.w[[1]]) %*% t(V.th[[1]])))
sigma.FFFI[[2]] <- V.th[[2]] %*% (Sigma.w[[2]] + sigma.b[2] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[2]]) + diag(diag(V.th[[2]] %*% (Sigma.v[[2]] - Sigma.w[[2]]) %*% t(V.th[[2]])))
sigma.FFFI[[3]] <- V.th[[3]] %*% (Sigma.w[[3]] + sigma.b[3] * (matrix(1, param[4], param[4]) - diag(param[4]))) %*% t(V.th[[3]]) + diag(diag(V.th[[3]] %*% (Sigma.v[[3]] - Sigma.w[[3]]) %*% t(V.th[[3]])))

# Well-separated clusters
mu.th.well <- rbind(rep(0, param[2]), rep(6, param[2]), rep(12, param[2]))
X.well <- vector(mode = "list", length = param[8])
for (s in 1:param[8]) {
  X.well[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.well[1, ], sigma.FFFI[[1]]),
                              mvrnorm(param[1]*param[6], mu.th.well[2, ], sigma.FFFI[[2]]),
                              mvrnorm(param[1]*param[7], mu.th.well[3, ], sigma.FFFI[[3]])),
                        nrow = param[1], ncol = param[2])
}

pugmm.FFFI.well <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FFFI.well <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.well[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FFFI.well <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FFFI.well[s, 1] <- adjustedRandIndex(label.th, pugmm.FFFI.well[[s]]$res$label)
  if (pugmm.FFFI.well[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FFFI.well[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FFFI.well[[s]]$res$model.name == "FFFI") {
    count.model <- count.model + 1
  }
  stat.FFFI.well[s, 2] <- pugmm.FFFI.well[[s]]$runtime[3]
}
mstat.FFFI.well <- data.frame(count.G, mean(stat.FFFI.well[, 1]), count.m, count.model, mean(stat.FFFI.well[, 2]))
colnames(mstat.FFFI.well) <- c("G", "mARI", "m", "model", "time")

# Close clusters
mu.th.close <- rbind(rep(0, param[2]), rep(3, param[2]), rep(7, param[2]))
X.close <- vector(mode = "list", length = 100)
for (s in 1:param[8]) {
  X.close[[s]] <- matrix(rbind(mvrnorm(param[1]*param[5], mu.th.close[1, ], sigma.FFFI[[1]]),
                               mvrnorm(param[1]*param[6], mu.th.close[2, ], sigma.FFFI[[2]]),
                               mvrnorm(param[1]*param[7], mu.th.close[3, ], sigma.FFFI[[3]])),
                         nrow = param[1], ncol = param[2])
}

pugmm.FFFI.close <- vector(mode = "list", length = param[8])
ncores <- detectCores()
cl <- makeCluster(ncores - 2)
registerDoParallel(cl)
system.time(
  pugmm.FFFI.close <- foreach (i = 1:param[8], .packages = c("PUGMM")) 
  %dopar% {
    s.time <- system.time( {
      set.seed(i)
      pugmm.res <- pugmm(scale(X.close[[i]]), G = 1:5, m = 1:5)
    })
    list(res = pugmm.res, runtime = s.time)
  }
)
stopCluster(cl)

stat.FFFI.close <- matrix(0, nrow = param[8], ncol = 2)
count.G <- 0
count.m <- 0
count.model <- 0
for (s in 1:param[8]) {
  stat.FFFI.close[s, 1] <- adjustedRandIndex(label.th, pugmm.FFFI.close[[s]]$res$label)
  if (pugmm.FFFI.close[[s]]$res$G == param[3]) {
    count.G <- count.G + 1
  }
  if (pugmm.FFFI.close[[s]]$res$m == param[4]) {
    count.m <- count.m + 1
  }
  if (pugmm.FFFI.close[[s]]$res$model.name == "FFFI") {
    count.model <- count.model + 1
  }
  stat.FFFI.close[s, 2] <- pugmm.FFFI.close[[s]]$runtime[3]
}
mstat.FFFI.close <- data.frame(count.G, mean(stat.FFFI.close[, 1]), count.m, count.model, mean(stat.FFFI.close[, 2]))
colnames(mstat.FFFI.close) <- c("G", "mARI", "m", "model", "time")

