This package provides a function for visualizing "response surfaces" of 
random forests as described in [Kumbier et al. 2018](https://arxiv.org/abs/#' 1810.07287). 
Briefly, response surfaces report the expected value of 
responses conditional on feature interactions identified by siRF. They are 
akin to partial dependency plots, with several key distinctions:
 
- We visualze expected response value on hold-out data rather than model 
predicted values.

- We extract targeted rules from the random forest model (matching signed 
interactions) and visualize response values as a function of these rules 
rather than over the full model.

## Installation

This package can be installed directly from github using devtools:

```
library(devtools)
devtools::install_github("karlkumbier/rfsurface")
```

`rfsurface` operates of fitted iRF models. The `iRF`` package can be 
similarly installed from github:

```
library(devtools)
devtools::install_github("karlkumbier/iRF")
```

[The vignette](vignette/vignette.html) provides a simple example.

