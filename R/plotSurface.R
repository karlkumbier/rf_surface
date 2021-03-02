#' Plot interaction
#'
#' Generate response surface plots for a given interaction.
#' @param x numeric feature matrix, with replicate features grouped
#' @param y response vector.  
#' @param int signed interaction to plot. Formatted as 'X1+_X2+_X3-_...'
#' @param varnames character vector indicating feature names. If NULL,
#'  colnames(x) are used as feature names.
#' @param read.forest output of readForest
#' @param qcut: quantile to define low/high levels of additional features beyond
#'  order-2 interations. Thresholds will generated using the specified quantile
#'  of random forest thresholds for corresponding features. 
#' @param col.pal color palette for response surfaces
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param zlab z-axis label
#' @param slab label for splitting variable
#' @param range.col range of response values for color palette
#' @param z.range z-axiz range
#' @param grid.surface size of grid to generate response surfaces over.
#' @param min.surface minimum number of observations required to generate a
#'  response surface.
#' @param min.nd minimum leaf node size to extract decision rules from.
#' @param pred.prob: if TRUE, z-axis indicates predicted probability from the
#'  random forest. If false, z-axis indicates distribution of responses y
#' @param plot.enrich: used for classification to plot response surface relative
#' to proportion of class-1 observation.
#' @param drop0: if TRUE, class-0 leaf nodes are removed before plotting
#' response surfaces.
#' @param nc: number of columns in plot
#' @param main plot title for response surfaces
#'
#' @export
#'
#' @importFrom rgl open3d persp3d par3d rgl.viewpoint movie3d spin3d mfrow3d 
#'  title3d
#' @importFrom dplyr select group_by summarize filter
#' @importFrom data.table data.table
#' @importFrom stringr str_split str_remove_all str_replace_all
#' @importFrom RColorBrewer brewer.pal
plotInt <- function(x, y, int, read.forest, 
                    varnames=NULL,
                    qcut=0.5,
                    col.pal=c('#1c3f66', '#306aab', 
                      '#6e96c4', '#ffb003',  '#ff8300'), 
                    xlab=NULL, ylab=NULL, zlab=NULL, slab=NULL,
                    range.col=NULL,
                    z.range=NULL,
                    grid.size=50,
                    min.surface=20,
                    min.nd=5,
                    pred.prob=FALSE,
                    plot.enrich=FALSE,
                    filt.rule=TRUE,
                    drop0=FALSE,
                    nc=2,
                    main=NULL) {
  
  # Check whether rgl package is installed
  if (! 'rgl' %in% rownames(installed.packages()))
    stop('Surface map plots require rgl installation')

  # Check whether read.forest is valid
  if (is.null(read.forest$node.feature))
    stop('read.forest missing node.feature')
  if (is.null(read.forest$node.obs))
    stop('read.forest missing node.obs')

  # Set feature names and check for replicates
  varnames <- groupVars(varnames, x)
  if (is.null(colnames(x))) {
    colnames(x) <- paste0('X', 1:ncol(x))
    varnames <- colnames(x)
  }
  
  if (any(duplicated(varnames))) 
    stop('Replicate features not supported')
  
  # Set z-axis scaling
  if (is.factor(y)) {    
    y <- as.numeric(y) - 1
  } 
  
  if (is.null(z.range)) {
    z.range <- range(y)
  }
  
  # Get feature indices for interaction in nf, x, and leaf ndoes
  int <- str_split(int, '_')[[1]]
  int.clean <- str_remove_all(int, '[-\\+]')
  int.nf <- int2Id(int.clean, varnames, split=TRUE, signed=FALSE)
  int.x <- int.nf %% ncol(x) + ncol(x) * (int.nf %% ncol(x) == 0)
  
  # Group node feature for interactions of any sign among given features
  p <- ncol(x)
  nf <- read.forest$node.feature[,1:p] + read.forest$node.feature[,(p + 1):(2 * p)]
  int.lf <- Matrix::rowMeans(nf[,int.nf] != 0) == 1

  if (length(int) > 2) {
    # Evaluate thresholds for features beyond order-2
    qCut <- function(x) quantile(x, probs = qcut)
    thr <- apply(nf[int.lf, int.nf], MAR=2, qCut)
    
    # Assign each observation to a unique plot 
    id.plot <- 3:length(int)
    id <- apply(t(x[,int.x[id.plot]]) > thr[id.plot], MAR=2, paste, collapse='_')
  } else {
    id <- rep(1, nrow(x))
    id.plot <- 1:2
  }
  
  # Generate grid of x/y values for surface maps
  grids <- quantileGrid(x, grid.size, int.x[1:2])
  
  # Extract hyperrectangles from RF decision paths
  rectangles <- forestHR(read.forest, int.nf, min.nd, int.lf)
  
  # Generate surface maps for each group of observations
  ids <- lapply(sort(unique(id)), '==', id)
  surf.scale <- ifelse(plot.enrich, mean(y), 1)

  surfaces <- lapply(ids, function(ii) {
    
    if (sum(ii) < min.surface) {
      warning('Fewer than min.surface observations, skipping surface')
      return(NULL)
    }

    genSurface(x[ii,], y[ii], int.nf[1:2], varnames=varnames, 
               rectangles=rectangles, min.nd=min.nd, surf.scale=surf.scale,
               filt.rule=filt.rule, drop0=drop0, grids=grids)
  })
  
  # Set color range for surface map plots
  if (is.null(range.col)) range.col <- range(unlist(surfaces))

  ngroup <- length(unique(id))
  if (ngroup > 4) 
    stop('Surface plots supported for up to order-4 interactions')
  if (ngroup == 1) {
    open3d()
    par3d(windowRect = c(0, 0, 1500, 1500))
  } 
  if (ngroup == 2) {
    open3d()
    nr <- ngroup / nc
    mfrow3d(nr=nr, nc=nc, sharedMouse = T)
    par3d(windowRect = c(0, 0, 1500, 1500))
  } 
  if (ngroup == 4) {
    open3d()
    nr <- ngroup / nc
    mfrow3d(nr=nr, nc=nc, sharedMouse = T)
    par3d(windowRect = c(0, 0, 1500, 1500))
  } 

  # Iterate over observation groups to generate response surfaces
  for (i in 1:length(unique(id))) {
    if (is.null(surfaces[[i]])) next
    
    # Generate title for surface map
    ii <- sort(unique(id))[i]
    ii <- str_replace_all(ii, 'TRUE', 'High')
    ii <- str_replace_all(ii, 'FALSE', 'Low')
    if (length(int) > 2) {
      i.split <- str_split(ii, '_')[[1]]
      if (is.null(slab)) slab <- int.clean[3:length(int)]
      xii <- paste(slab, i.split, sep=': ')
      main.ii <- paste(main, paste(xii, collapse=', '), collapse=' - ')
    } else {
      main.ii <- main
    } 

    # Generate response surface for curent group
    plotInt2(surfaces[[i]], xlab=xlab, ylab=ylab, zlab=zlab, main=main.ii,
             col.pal=col.pal, range.col=range.col, z.range=z.range)
    rgl.viewpoint(zoom=0.95, theta=-5, phi=-60)
  
  }
}

plotInt2 <- function(surface,
                     col.pal=c('#1c3f66', '#306aab', 
                       '#6e96c4', '#ffb003',  '#ff8300'), 
                     xlab=NULL, ylab=NULL, zlab=NULL,
                     main=NULL,
                     range.col=NULL,
                     filt.rule=TRUE,
                     z.range=c(0,1),
                     n.cols=100, 
                     axes=TRUE) {
  # Generates surface map plot of order-2 interaction
  # args:
  #   surface: response surface matrix, output of genSurface
  #   col.pal: color palette of surface map
  #   xlab, ylab, zlab: axis labels
  #   main: title for the plot
  #   range.col: range of color values
  #   z.range: range for response axis
  #   n.cols: number of colors in color pal
  #   axes: T/F indicating whether axes should be plotted
  range.col <- c(range.col, as.vector(surface))
  if (length(unique(range.col)) == 1) n.cols <- 1
  
  palette <- colorRampPalette(col.pal)
  if (length(unique(range.col)) == 1) {
    colors <- palette(1)
    facet.col <- 1
  } else {
    colors <- palette(n.cols)
    facet.col <- cut(range.col, n.cols)[-seq(2)]
  }

  # Plot interaction response surface
  par3d(cex=1.5)
  if (is.null(zlab)) zlab <- ''
  g1n <- as.numeric(rownames(surface))
  g2n <- as.numeric(colnames(surface))
  persp3d(x=g1n, y=g2n, z=surface, xlab=xlab, ylab=ylab, zlab=zlab, 
          zlim=z.range, col=colors[facet.col], axes=axes)
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
