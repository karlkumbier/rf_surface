#' Plot interaction
#'
#' Visualization of response surface plots.
#' @param x numeric feature matrix, with replicate features grouped
#' @param int signed interaction to plot. If numeric, int is assumed to
#'   correspond to column indices to be plotted for interaction. If character,
#'   assumed to be formatted as 'X1+_X2+_X3-_...'
#' @param y response vector.  
#' @param fit a fitted random forest, from packages randomForest or ranger.
#' @param read.forest output of readForest
#' @param varnames character vector indicating feature names. If NULL,
#'  colnames(x) are used as feature names.
#' @param col.pal color palette for response surfaces. A function that takes an
#'   integer input and returns colors to be use in the palette.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param zlab z-axis label
#' @param slab for order 3 and 4 interactions, label for split plots
#' @param z.range z-axiz range
#' @param nbin: number of bins to plot surface map over
#' @param min.surface minimum number of observations required to generate a
#'  response surface.
#' @param filter.rules: a list of filtering functions to be applied to rf
#'   decision paths. If NULL, default rules will filter to a random sample of
#'   10% of leaf nodes with at least 5 observations.
#' @param filter_x a filtering function to be applied to data matrix. Takes as
#'   arguments x (data matrix), int (numeric vector of interaciton ids), and
#'   thresholds (numeric vector of rf thresholds, columns corresponding to
#'   features in int), returns indices of x to be kept.
#' @param wt.node indicator for how nodes are to be weighted in response
#'   surfaces. One of `size` - weighting proportional to leaf node size or
#'   `none` - indicating uniform weighting.
#' @param type one of `rgl` - 3d response surface or ggplot - 2d response
#'   surface
#' @param main plot title for response surfaces
#'
#' @export
#'
#' @importFrom rgl open3d persp3d par3d rgl.viewpoint movie3d spin3d title3d
#' @importFrom dplyr select group_by summarize filter
#' @importFrom data.table data.table
#' @importFrom stringr str_split str_remove_all str_replace_all
#' @importFrom viridis magma
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradientn xlab ylab labs
plotInt <- function(x, int,  
                     y=NULL,
                     fit=NULL,
                     read.forest=NULL,
                     varnames=colnames(x),
                     col.pal=magma, 
                     xlab=NULL, 
                     ylab=NULL, 
                     zlab=NULL, 
                     slab=NULL,
                     z.range=NULL,
                     nbin=50,
                     filter.rules=NULL,
                     filterX=NULL,
                     wt.node='size',
                     type='rgl',
                     main=NULL) {
 
  n <- nrow(x)
  p <- ncol(x)
  pred.prob <- is.null(y)

  # Check for one of read.forest/fit
  if (is.null(read.forest) & is.null(fit)) {
      stop('Specify one of `read.forest` or `fit`')
  }

  # Read out RF decision paths
  if (is.null(read.forest)) {
    read.forest <- readForest(fit, x=x, oob.importance=FALSE)
  }

  # Check whether read.forest is valid
  if (is.null(read.forest$node.feature)) stop('read.forest missing node.feature')
  if (is.null(read.forest$node.obs)) stop('read.forest missing node.obs')
  
  # Set feature names and check for replicates
  varnames <- groupVars(varnames, x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0('X', 1:ncol(x))
    varnames <- colnames(x)
  }
  
  # Check for duplicate features
  if (any(duplicated(varnames))) stop('Replicate features not supported')

  # Convert binary factor
  if (is.factor(y)) y <- as.numeric(y) - 1
  
  # Set z-axis scaling
  if (is.null(z.range) & !pred.prob) z.range <- range(y)
  if (is.null(z.range) & pred.prob) z.range <- range(read.forest$tree.info$prediction)
 
  # Check for valid interaction and convert to numeric IDs
  if (!is.numeric(int)) {
      signed <- str_detect(int, '(\\+|-)')
      int <- int2Id(int, varnames, signed=signed)
      int <- int %% p + p * (int %% p == 0)
  }
  
  # Collapse node feature matrix to unsigned
  if (ncol(read.forest$node.feature) == 2 * p) {
    read.forest$node.feature <- read.forest$node.feature[,1:p] + 
      read.forest$node.feature[,(p + 1):(2 * p)]
  }
  
  # Generate grid of x/y values for surface maps
  bins <- quantileGrid(x, nbin, int[1:2])
  
  # Extract hyperrectangles from RF decision paths
  rectangles <- forestHR(read.forest, int)
  
  # Filter data matrix if rules specified
  if (!is.null(filterX)) {
    id <- filterX(x, int, rectangles$splits)
    x <- x[id,]
    if (!is.null(y)) y <- y[id]
  }

  # Generate surface for current plot
  surface <- genSurface(x, int[1:2],
                        y=y,
                        varnames=varnames, 
                        rectangles=rectangles, 
                        wt.node=wt.node,
                        filter.rules=filter.rules,
                        bins=bins)
  
  # Set quantile names for grid
  colnames(surface) <- seq(0, 1, length.out=nrow(surface))
  rownames(surface) <- seq(0, 1, length.out=ncol(surface))
  
  # Select plotting method, one of rgl or ggplot
  plotFun <- ifelse(type == 'rgl', rglplotSurface2, ggplotSurface2)
  
  # Generate response surface for curent group
  plotFun(surface, 
          col.pal=col.pal,
          xlab=xlab, 
          ylab=ylab, 
          zlab=zlab, 
          main=main,
          z.range=z.range)
  
}


ggplotSurface2 <- function(surface,
                           col.pal=magma, 
                           xlab=NULL, 
                           ylab=NULL, 
                           zlab=NULL,
                           main=NULL,
                           z.range=range(surface)) {

  # TODO: fixed z-range for visualizations
  # Set axis names
  xlab <- ifelse(is.null(xlab), '', xlab)
  ylab <- ifelse(is.null(ylab), '', ylab)
  zlab <- ifelse(is.null(zlab), '', zlab)
  
  p <- reshape2::melt(surface) %>%
    ggplot(aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradientn(colours=col.pal(100)) +
    xlab(xlab) +
    ylab(ylab) +
    labs(fill=zlab)
  
  plot(p)
}

rglplotSurface2 <- function(surface,
                            col.pal=magma, 
                            xlab=NULL, 
                            ylab=NULL, 
                            zlab=NULL,
                            main=NULL,
                            z.range=range(surface),
                            axes=TRUE) {
  # Generates surface map plot of order-2 interaction
  # args:
  #   surface: response surface matrix, output of genSurface
  #   col.pal: color palette of surface map
  #   xlab, ylab, zlab: axis labels
  #   main: title for the plot
  #   z.range: range for response axis
  #   axes: T/F indicating whether axes should be plotted
  
  # Initialize color palette
  if (length(unique(c(surface))) == 1) {
    facet.col <- 1
  } else {
    facet.col <- as.numeric(cut(c(surface), 100))
  }
  
  # Generate color palette for response surface
  #n.cols <- min(100, length(unique(facet.col)))
  colors <- col.pal(100)
  
  # Set axis names
  xlab <- ifelse(is.null(xlab), '', xlab)
  ylab <- ifelse(is.null(ylab), '', ylab)
  zlab <- ifelse(is.null(zlab), '', zlab)
  
  # Plot interaction response surface
  par3d(cex=1.5)
  persp3d(x=as.numeric(rownames(surface)), 
          y=as.numeric(colnames(surface)), 
          z=surface, 
          xlab=xlab, 
          ylab=ylab, 
          zlab=zlab, 
          zlim=z.range, 
          col=colors[facet.col], 
          axes=axes)
  
  if (!is.null(main)) title3d(main, line = 3)
}

quantileGrid <- function(x, nbin, int) {
  # Generate a grid at quantiles of x
  stopifnot(length(int) == 2)
  nbin <- nbin + 1
  bins <- list()
  bins$g1 <- quantile(x[,int[1]], probs=seq(0, 1, length.out=nbin))
  bins$g2 <- quantile(x[,int[2]], probs=seq(0, 1, length.out=nbin))
  return(bins)
}
