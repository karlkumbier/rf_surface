library(iRF)
library(tidyverse)
library(viridis)

source('~/github/Rpackages/rf_surface/R/utilities.R')
source('~/github/Rpackages/rf_surface/R/generateSurface.R')
source('~/github/Rpackages/rf_surface/R/plotSurface.R')

set.seed(47)
n <- 1000
p <- 50

x <- matrix(rnorm(n * p), nrow=n, ncol=p)
colnames(x) <- str_c('X', 1:ncol(x))
y <- as.numeric(x[,1] > 0 & x[,2] > 0)

# Random swapping noise
id.swap <- sample(n, 100)
y[id.swap] <- 1 - y[id.swap]

if (!file.exists('fit.Rdata')) {
  fit <- iRF(x=x, y=as.factor(y), n.iter=2, type='ranger')
  read.forest <- readForest(fit$rf.list, x=x, oob.importance=FALSE)
  save(file='fit.Rdata', fit, read.forest)
} else {
  load('fit.Rdata')
}

# Generate response surface for select rule list
rule.list <- list()
rule.list[[1]] <- function(x) sample_frac(x, 0.1)

plotInt(x=x, int='X1+_X2+', read.forest=read.forest, y=y, z.range=0:1)

