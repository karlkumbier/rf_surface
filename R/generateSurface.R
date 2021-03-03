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
#TODO: return..., imports - subsetReadForest, ??? filt.rule, surf.scale, drop0 ???
#TODO: is int assumed indices or interaction format?
#TODO: update largest.rule to be one of largest, sample, full (and class 0/1),
#   general filtering function
genSurface <- function(x, y, int, 
                       read.forest=NULL,
                       rectangles=NULL,
                       varnames=NULL,
                       nbin=100,
                       bins=NULL,
                       pred.prob=FALSE,
                       largest.rule=TRUE,
                       surf.scale=1,
                       drop0=FALSE,
                       min.nd=5) {
  
  n <- nrow(x)
  p <- ncol(x)
 
  # Check for valid interaction
  if (!is.numeric(int)) {
      # TODO: check that this doesn't create issues
      # TODO: support for unsigned interactions?
      int <- int2Id(int, varnames, signed=TRUE)
      int <- int %% p + p * (int == p)
  }
  
  if (length(int) != 2)
    stop('Response surface can only be generated over 2 features')

  # Check for one of read.forest/rectangles
  if (is.null(read.forest) & is.null(rectangles))
      stop('Specify one of `rectangles` or `read.forest`')

  # Set feature names and check for replicates
  varnames <- groupVars(varnames, x)
  if (any(duplicated(varnames)))
    stop('Replicate features not supported')

  # Extract hyperrectangles from readForest output
  if (is.null(rectangles)) { 
      # Convert standard format int to indices
      rectangles <- forestHR(read.forest, int, min.nd)
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
  rules <- list()
  rules[[1]] <- function(x) filter(x, prediction == 1)
  rules[[2]] <- function(x) group_by(x, tree) %>% filter(size.node == max(size.node))

  # Filter leaf node hyperrectangles
  rectangles <- filterHR(rectangles, rules)

  # Get thresholds and sign for interaction rules
  tt <- rectangles$splits

  # TODO: RCPP wrapper for this
  # Evaluate distriution of responses across each decision rule
  grid <- matrix(0, nrow=nbin, ncol=nbin)

  if (wt.node == 'none') wt <- rep(1, nrow(tt))
  if (wt.node == 'size.node') wt <- rectangles$nodes$size.node

  for (i in 1:nrow(tt)) {
    
    # weight response surfaces based on size of leaf node

    # Evalaute which observations/grid elements correspond to current HR
    idcs1 <- g1 >= tt[i, 1]
    x1 <- x[,int[1]] >= tt[i, 1]

    idcs2 <- g2 >= tt[i, 2]
    x2 <- x[,int[2]] >= tt[i, 2]

    if (pred.prob) {
      # Evaluate RF predictions for region corresponding to current HR
      yy <- rectangles$prediction[i]
      grid[idcs1, idcs2] <- grid[idcs1, idcs2] + yy * wt[i]
    } else {
      if (any(x1 & x2))
        grid[idcs1, idcs2] <- grid[idcs1, idcs2] +  mean(y[x1 & x2]) * wt[i]
      if (any(!x1 & x2))
        grid[!idcs1, idcs2] <- grid[!idcs1, idcs2] +  mean(y[!x1 & x2]) * wt[i]
      if (any(x1 & !x2))
        grid[idcs1, !idcs2] <- grid[idcs1, !idcs2] +  mean(y[x1 & !x2]) * wt[i]
      if (any(!x1 & !x2))
        grid[!idcs1, !idcs2] <- grid[!idcs1, !idcs2] +  mean(y[!x1 & !x2]) * wt[i]
    }
  }

  # Rescale surface for node size and generate corresponding color palette
  nsurface <- sum(rectangles$size.node)
  if (nsurface != 0) grid <- grid / nsurface
  if (all(grid == 0)) grid <- grid + 1e-3
  rownames(grid) <- g1n
  colnames(grid) <- g2n
  return(grid)
}

filterHR <- function(rectangles, rules) {
    # Applies a collection of filtering functions to leaf nodes
    if (!is.list(rules)) rukes <- list(rules)
    
    # Iterate over filter functions
    for (r in rules) {
       rectangles$nodes <- r(rectangles$nodes) 
    }

    # Subset threshold matrix based on remaining nodes
    rectangles$splits <- rectangles$splits[rectangles$nodes$ID,]
    return(rectangles)
}

forestHR <- function(read.forest, int, min.nd) {
  # Read hyperrectangles from RF for a specified interactin
  # args:
  #   read.forest: list as returned by readForest, including node.feature 
  #     and tree.info entries
  #   int: vector of indices specifying features for hyperrectangles
  #   min.nd: minimum node size to extract hyperrectangles from
  
  # Extract splitting features and thresholds
  p <- ncol(read.forest$node.feature) / 2
  read.forest$node.feature <- read.forest$node.feature[,1:p] + read.forest$node.feature[,(p + 1):(2 * p)]
 
  # Set active nodes for interaction
  int.lf <- Matrix::rowMeans(read.forest$node.feature[,int] != 0) == 1

  # Subset to active nodes larger than min.nd
  id.subset <- int.lf & read.forest$tree.info$size.node >= min.nd
  read.forest <- subsetReadForest(read.forest, id.subset)

  
  # Group data by leaf node for return
  nodes <- read.forest$tree.info %>% 
      select(prediction, node.idx, tree, size.node) %>%
      mutate(ID=1:n())

  splits <- read.forest$node.feature[,int]

  return(list(nodes=nodes, splits=splits, int=int))
}
