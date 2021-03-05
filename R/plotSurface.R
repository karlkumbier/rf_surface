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
plotInt <- function(x, int, read.forest,
                    y=NULL,
                    varnames=colnames(x),
                    quantile.cut=0.5,
                    col.pal=magma, 
                    xlab=NULL, 
                    ylab=NULL, 
                    zlab=NULL, 
                    slab=NULL,
                    z.range=NULL,
                    nbin=50,
                    min.surface=20,
                    min.nd=5,
                    pred.prob=FALSE,
                    filter.rules=NULL,
                    wt.node='size',
                    type='rgl',
                    main=NULL) {
  
  n <- nrow(x)
  p <- ncol(x)
  
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
  if (is.null(z.range) & !is.null(y)) z.range <- range(y)
  if (is.null(z.range) & is.null(y)) 
    z.range <- range(read.forest$tree.info$prediction)
  
  # Get feature indices for interaction in nf, x, and leaf ndoes
  int <- str_remove_all(str_split(int, '_')[[1]], '[-\\+]')
  int <- int2Id(int, varnames, split=TRUE, signed=FALSE)
  
  # Collapse node feature matrix to unsigned
  nf <- read.forest$node.feature
  if (ncol(nf) == 2 * p) nf <- nf[,1:p] + nf[,(p + 1):(2 * p)]
  
  # Determine leaf nodes containing interaction
  int.lf <- Matrix::rowMeans(nf[,int] != 0) == 1
  
  if (length(int) > 2) {
    # Evaluate thresholds for features beyond order-2
    qcut <- function(x) quantile(x, probs = quantile.cut)
    thr <- apply(nf[int.lf, int], MAR=2, qcut)
    
    # Assign each observation to a unique plot 
    id.plot <- 3:length(int)
    x.thresh <- t(x[,int[id.plot]]) > thr[id.plot]
    id <- apply(x.thresh, MAR=2, paste, collapse='_')
  } else {
    id <- rep(1, nrow(x))
    id.plot <- 1:2
  }
  
  # Generate grid of x/y values for surface maps
  bins <- quantileGrid(x, nbin, int[1:2])
  
  # Extract hyperrectangles from RF decision paths
  rectangles <- forestHR(read.forest, int)
  
  # Generate surface maps for each group of observations
  surfaces <- lapply(sort(unique(id)), function(iid) {
    
    # Get obserations corresponding to current plot
    ii <- id == iid
    
    # Check that the current plot contains the minimum number of observations
    if (sum(ii) < min.surface) {
      warning('Fewer than min.surface observations, skipping surface')
      return(NULL)
    }
    
    # Set response to be passed in for surface
    if (is.null(y)) {
      ysurface.i <- NULL
    } else {
      ysurface.i <- y[ii]
    }
    
    # Generate surface for current plot
    genSurface(x[ii,], int[1:2],
               y=ysurface.i,
               varnames=varnames, 
               rectangles=rectangles, 
               wt.node=wt.node,
               pred.prob=pred.prob,
               filter.rules=filter.rules,
               bins=bins)
  })
  
  # Initialize windows for rgl plots
  if (type == 'rgl') {
    ngroup <- length(unique(id))
    if (ngroup > 4) 
      stop('Surface plots supported for up to order-4 interactions')
    
    # Initialize window for 1 response surface, order-2 interaction
    if (ngroup == 1) {
      open3d()
      par3d(windowRect = c(0, 0, 1500, 1500))
    } 
    
    # Initialize window for 2 response surface, order-3 interaction
    if (ngroup == 2) {
      open3d()
      nr <- ngroup / 2
      mfrow3d(nr=nr, nc=2, sharedMouse = T)
      par3d(windowRect = c(0, 0, 1500, 1500))
    } 
    
    # Initialize window for 2 response surface, order-4 interaction
    if (ngroup == 4) {
      open3d()
      nr <- ngroup / 2
      mfrow3d(nr=nr, nc=2, sharedMouse = T)
      par3d(windowRect = c(0, 0, 1500, 1500))
    } 
  }
  
  #TODO: quantile names for grid
  
  # Iterate over observation groups to generate response surfaces
  for (i in 1:length(unique(id))) {
    if (is.null(surfaces[[i]])) next
    
    # Set quantile names for grid
    colnames(surfaces[[i]]) <- seq(0, 1, length.out=nbin + 1)
    rownames(surfaces[[i]]) <- seq(0, 1, length.out=nbin + 1)
    
    
    
    # Reformat group names for title
    ii <- sort(unique(id))[i]
    ii <- str_replace_all(ii, 'TRUE', 'High')
    ii <- str_replace_all(ii, 'FALSE', 'Low')
    
    # Generate title for response surface
    if (length(int) > 2) {
      i.split <- str_split(ii, '_')[[1]]
      if (is.null(slab)) slab <- int.clean[3:length(int)]
      xii <- paste(slab, i.split, sep=': ')
      main.ii <- paste(main, paste(xii, collapse=', '), collapse=' - ')
    } else {
      main.ii <- main
    } 
    
    # Select plotting method, one of rgl or ggplot
    plotFun <- ifelse(type == 'rgl', rglplotInt2, ggplotInt2)
    
    # Generate response surface for curent group
    plotFun(surfaces[[i]], 
            col.pal=col.pal,
            xlab=xlab, 
            ylab=ylab, 
            zlab=zlab, 
            main=main.ii)
    
    if (type == 'rgl') rgl.viewpoint(zoom=0.95, theta=-5, phi=-60)
  }
}


ggplotInt2 <- function(surface,
                       col.pal=magma, 
                       xlab=NULL, 
                       ylab=NULL, 
                       zlab=NULL,
                       main=NULL,
                       z.range=range(surface),
                       axes=NULL) {
  
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

rglplotInt2 <- function(surface,
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
