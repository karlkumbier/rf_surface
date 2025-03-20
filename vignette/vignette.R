#' # Overview
#' This package provides a function for visualizing "response surfaces" of 
#' random forests as described in [Kumbier et al. 2018](https://arxiv.org/abs/#' 1810.07287). Briefly, response surfaces report the expected value of 
#' responses conditional on feature interactions identified by siRF. They are 
#' akin to partial dependency plots, with several key distinctions:
#'  
#' - We visualze expected response value on hold-out data rather than model 
#' predicted values.
#' 
#' - We extract targeted rules from the random forest model (matching signed 
#' interactions) and visualize response values as a function of these rules 
#' rather than over the full model.
#' 
#' ## Installation
#' 
#' This package can be installed directly from github using devtools:
#' 
#' ```
#' library(devtools)
#' devtools::install_github("karlkumbier/rfsurface")
#' ```
#' 
#' `rfsurface` operates of fitted iRF models. The `iRF`` package can be 
#' similarly installed from github:
#' 
#' ```
#' library(devtools)
#' devtools::install_github("karlkumbier/iRF")
#' ``` 
#'
#' # Example
#' The example below demonstrates how to visualize response surfaces from a 
#' fitted iRF. Responses are generated as random Bernoulli variables based on #' the decision rule: $X_1 > median(X_1) \& X_2 > median(X_2)$. The sampling 
#' probability for responses is set to 0.8 when the rule is active and 0.2 
#' otherwise. We include an additional 48 features in the dataset as noise.
#+ example, message=FALSE, warning=FALSE
library(rfsurface)
library(iRF)
library(tidyverse)

set.seed(47)
n <- 2000
p <- 50

x <- data.frame(matrix(rexp(n * p), nrow=n, ncol=p))
yrule <- as.numeric(x[,1] > median(x[,1]) & x[,2] > median(x[,2]))
yprob <- 0.8 * yrule + 0.2 * (1 - yrule)
y <- rbinom(n, 1, yprob)

# Set training and test indices
train.idx <- sample(1:n, size=0.5*n)
test.idx <- setdiff(1:n, train.idx)

# Fit iRF model
fit <- iRF(
    x=x[train.idx,], 
    y=as.factor(y)[train.idx], 
    n.iter=2, 
    type='ranger',
    int.return=2
)

interactions <- filter(fit$interaction, stability > 0.5)
print(head(interactions))

#' iRF identifies the interaction `X_1+_X_2+` as both stable and predictive. 
#' Below, we extract decision rules from the fitted model and visualize test 
#' set response distribution as a function of these rules. Specifically, each 
#' rule of the form $X_1 > t_1 \& X_2 > t_2$ corresponds to a hyperrectangle. 
#' We plot the expected value of responses within and outside of this 
#' hyperrectangle.
#' 
#' Note that in our simulation, $X_i$ were drawn as random exponential 
#' variables. We can control the x and y-axis scaling using the argument 
#' `binFun`, which we set to $log(x)$ below.
#+ plot, message=FALSE, warning=FALSE
# Extract rules from fitted iRF model
read.forest <- readForest(fit$rf.list, x=x[test.idx,], oob.importance=FALSE)

# Generate response surface for select rule list
plotInt(
    x=x[test.idx,], 
    int='X1+_X2+', 
    read.forest=read.forest, 
    y=y[test.idx], 
    z.range=0:1, 
    binFun=function(x) log(x)
)

