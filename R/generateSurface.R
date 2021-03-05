#' Generate response surface for 2D interaction
#'
#' Response surface give E(Y|X) over regions learned by an iRF.
#'
#' @param x numeric feature matrix, with replicate features grouped
#' @param y response vector.  
#' @param int signed interaction to plot. Formatted as 'X1+_X2+_X3-_...'
#' @param read.forest output of readForest.
#' @param rectangles a list of hyperrectangles corresponding to leaf nodes in an
#'  RF, as retuned by forestHR. If both rectangles and read.forest are supplied,
#'  read.forest will be ignored.
#' @param varnames character vector indicating feature names. If NULL,
#'  colnames(x) are used as feature names.
#' @param nbin: number of bins to plot surface map over
#' @param bins: user generated grid to plot over. If supplied, nbin is 
#'  ignored
#' @param pred.prob: T/F indicating whether the z axis should indicate predicted
#'  probability from the fitted RF (T) or average value of y (F). 
#'  TODO: if y=NULL, default to T.
#' @param largest.rule: should response surfaces be generated using only the
#'  largest rule per tree? Plotting response surfaces over a subset of leaf
#'  nodes can substantially reduce computation time.
#' @param wt.node: one of 'none' or 'size' - indicating how leaf nodes are
#'  weightied in the response surface. none = equal weighting, size = leaf nodes
#'  weighted proportional to number of observations.
#' @param min.nd: minimum leaf node size to used in generating response surface.
#'  Leaf nodes containing fewer than min.nd observations in fitted forest will
#'  be dropped before generating the response surface.
#'
#' @export
genSurface <- function(x, int,
                       y=NULL,
                       read.forest=NULL,
                       rectangles=NULL,
                       wt.node='size',
                       varnames=colnames(x),
                       nbin=100,
                       bins=NULL,
                       pred.prob=FALSE,
                       filter.rules=NULL) {
  
  n <- nrow(x)
  p <- ncol(x)
 
  if (is.null(y) & !pred.prob) {
    pred.prob <- TRUE
    warning('No response supplied, using pred.prob=TRUE')
  }
  
  # Check for valid interaction
  if (!is.numeric(int)) {
      signed <- str_detect(int, '(\\+|-)')
      int <- int2Id(int, varnames, signed=signed)
      int <- int %% p + p * (int == p)
  }
  
  if (length(int) != 2) {
    stop('Response surface can only be generated over 2 features')
  }

  # Check for one of read.forest/rectangles
  if (is.null(read.forest) & is.null(rectangles)) {
      stop('Specify one of `rectangles` or `read.forest`')
  }

  # Set feature names and check for replicates
  varnames <- groupVars(varnames, x)
  if (any(duplicated(varnames))) {
    stop('Replicate features not supported')
  }
  
  # Extract hyperrectangles from readForest output
  if (is.null(rectangles)) { 
    
    # Collapse node feature matrix - will use any int between features, 
    # regardless of sign
    if (ncol(read.forest$node.feature) == 2 * p) {
      read.forest$node.feature <- read.forest$node.feature[,1:p] + 
        read.forest$node.feature[,(p + 1):(2 * p)]
    }
    
    # Convert standard format int to indices
    rectangles <- forestHR(read.forest, int)
  }

  # Generate grid to plot surface over either as raw values or quantiles
  if (is.null(bins)) {
    # Sequences spanning the range of each interacting feature
    g1 <- seq(min(x[,int[1]]), max(x[,int[1]]), length.out=nbin)
    g2 <- seq(min(x[,int[2]]), max(x[,int[2]]), length.out=nbin)
  } else {
    # Pre-specified grid sequence
    g1 <- bins$g1
    g2 <- bins$g2
    nbin <- length(g1)
  }
  
  g1n <- round(g1, 2)
  g2n <- round(g2, 2)

  # Define a set of functions for filtering leaf node hyperrectangles
  if (is.null(filter.rules)) {
    filter.rules <- list()
    
    filter.rules[[1]] <- function(x) {
      filter(x, size.node >= 5) 
    }
    
    filter.rules[[2]] <- function(x) {
      group_by(x, tree) %>% sample_n(1) 
    }
  }

  # Filter leaf node hyperrectangles
  rectangles <- filterHR(rectangles, filter.rules)

  # Get thresholds and sign for interaction rules
  thresholds <- rectangles$splits

  # Evaluate distriution of responses across each decision rule
  grid <- matrix(0, nrow=nbin, ncol=nbin)

  if (wt.node == 'none') wt <- rep(1, nrow(thresholds))
  if (wt.node == 'size') wt <- rectangles$nodes$size.node

  nsurface <- 0
  for (i in 1:nrow(thresholds)) {
    
    # weight response surfaces based on size of leaf node

    # Evalaute which observations/grid elements correspond to current HR
    i1 <- g1 >= thresholds[i, 1]
    x1 <- x[,int[1]] >= thresholds[i, 1]

    i2 <- g2 >= thresholds[i, 2]
    x2 <- x[,int[2]] >= thresholds[i, 2]

    if (pred.prob) {
      # Evaluate RF predictions for region corresponding to current HR
      y <- rectangles$nodes$prediction[i]
      grid[i1, i2] <- grid[i1, i2] + y * wt[i]
    } else {
      
      if (any(x1 & x2)) {
        nsurface <- nsurface + wt[i]
        grid[i1, i2] <- grid[i1, i2] +  mean(y[x1 & x2]) * wt[i]
        
        y.inactive <- mean(y[!x1 | !x2]) * wt[i]
        grid[!i1, i2] <- grid[!i1, i2] +  y.inactive
        grid[i1, !i2] <- grid[i1, !i2] +  y.inactive
        grid[!i1, !i2] <- grid[!i1, !i2] +  y.inactive
      }

    }
  }

  # Rescale surface for node size
  if (nsurface != 0) grid <- grid / nsurface
  if (all(grid == 0)) grid <- grid + 1e-10
  
  rownames(grid) <- g1n
  colnames(grid) <- g2n
  return(grid)
}

filterHR <- function(rectangles, rules) {
    # Applies a collection of filtering functions to leaf nodes
    if (!is.list(rules)) rules <- list(rules)
    
    # Iterate over filter functions
    for (r in rules) {
       rectangles$nodes <- r(rectangles$nodes) 
    }

    # Subset threshold matrix based on remaining nodes
    rectangles$splits <- rectangles$splits[rectangles$nodes$ID,]
    return(rectangles)
}

forestHR <- function(read.forest, int) {
  # Read hyperrectangles from RF for a specified interactin
  # args:
  #   read.forest: list as returned by readForest, including node.feature 
  #     and tree.info entries
  #   int: vector of indices specifying features for hyperrectangles

  # Set active nodes for interaction
  int.lf <- Matrix::rowMeans(read.forest$node.feature[,int] != 0) == 1
  
  # Group data by leaf node for return
  nodes <- read.forest$tree.info %>% 
      select(prediction, node.idx, tree, size.node) %>%
      mutate(ID=1:n())

  splits <- read.forest$node.feature[,int]

  return(list(nodes=nodes, splits=splits, int=int))
}
