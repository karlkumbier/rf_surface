library(iRF)
library(tidyverse)
library(rgl)
library(viridis)
source('~/github/Rpackages/rf_surface/R/utilities.R')
source('~/github/Rpackages/rf_surface/R/generateSurface.R')
source('~/github/Rpackages/rf_surface/R/plotSurface.R')

n <- 1000
p <- 50

x <- matrix(rnorm(n * p), nrow=n, ncol=p)
colnames(x) <- str_c('X', 1:ncol(x))
y <- as.numeric(x[,1] > 0 & x[,2] > 0)

# Random swapping noise
id.swap <- sample(n, 100)
y[id.swap] <- 1 - y[id.swap]

fit <- iRF(x=x, y=as.factor(y), n.iter=2, type='ranger')
read.forest <- readForest(fit$rf.list, x=x, oob.importance=FALSE)

# Generate response surface for select rule list
rule.list <- list()
rule.list[[1]] <- function(x) group_by(x, tree) %>% sample_frac(0.1)
surface <- genSurface(x, 'X1+_X2+', y=y, read.forest=read.forest, filter.rules=rule.list)
plotInt(x, 'X1+_X2+', read.forest, y)
plotInt(x, 'X1+_X2+', read.forest, y, type='ggplot', xlab='x', ylab='y', zlab='p')
x

