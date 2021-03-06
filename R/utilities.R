#' @importFrom stringr str_replace_all str_remove_all str_split
groupVars <- function(varnames.grp, x) {
  # Generate grouped variable names
  if (is.null(varnames.grp)) {
    if (!is.null(colnames(x))) {
      varnames.grp <- colnames(x)
    } else {
      varnames.grp <- 1:ncol(x)
    }
  }

  stopifnot(length(varnames.grp) == ncol(x))
  return(varnames.grp)
}

pasteInt <- function(x) {
  # Combine interaction into single string
  x <- paste(x, collapse='_')
  return(x)
}

nameInts <- function(ints, varnames, signed=TRUE) {
  # Convert interactions indicated by indices to interactions indicated by
  # variable names. Naming convention for an interaction is:
  #   <variable1(sign)>_ <variable2(sign)>_...

  varnames <- unique(varnames)
  p <- length(varnames)
  if (signed)
    signs <- lapply(ints, function(z) ifelse(z > p, '+', '-'))
  else
    signs <- ''

  # Adjust indexing to match varnames
  ints <- lapply(ints, function(z) z %% p + p * (z == p | z == 2 * p))
  ints.name <- mapply(function(i, s) nameInt(varnames, i, s), ints, signs)
  return(ints.name)
}

nameInt <- function(varnames, idx, sgn) {
  int <- paste0(varnames[idx], sgn)
  int <- paste(sort(int), collapse='_')
  return(int)
}

int2Id <- function(int, varnames.grp, signed=FALSE, split=FALSE) {
  # Determine integer index of named variable (signed or not)
  if (!split) int <- str_split(int, '_')[[1]]

  if (signed) {
    sgn <- grep('\\+$', int)
    varnames.grp <- str_remove_all(varnames.grp, '[\\+\\-]')
    int <- str_remove_all(int, '[\\+\\-]')
  }

  varnames.grp <- unique(varnames.grp)
  id <- sapply(int, function(i) which(varnames.grp == i))
  if (signed) {
    adjust <- rep(0, length(int))
    adjust[sgn] <- length(varnames.grp)
    id <- id + adjust
  }

  return(id)
}

unsign <- function(int) {
  # Remove sign indicators from interaction strings
  return(str_replace_all(as.character(int), '[-\\+]', ''))
}

intSign <- function(int, split=TRUE) {
  # Evaluate sign of interactions
  if (!split) int <- str_split(int, '_')[[1]]
  sgn <- rep(-1, length(int))
  sgn[grep('\\+', int)] <- 1
  return(sgn)
}

intSubsets <- function(int, split=TRUE) {
  # Generate order 1, s - 1, and s subsets of an order-s interaction
  if (!split) int <- str_split(as.character(int), '_')[[1]]
  if (length(int) == 1) return(int)
  sub.ord <- c(1, length(int) - 1, length(int))
  subs <- lapply(sub.ord, combn, x=int, simplify=FALSE)
  subs <- unlist(subs, recursive=FALSE)
  return(subs)
}

lreplicate <- function(n, expr, ...) {
  # replicate with list return
  out <- replicate(n, expr, ..., simplify=FALSE)
  return(out)
}

subsetReadForest <- function(read.forest, subset.idcs) {
  # Subset nodes from readforest output
  if (!is.null(read.forest$node.feature))
    read.forest$node.feature <- read.forest$node.feature[subset.idcs,]

  if (!is.null(read.forest$tree.info))
    read.forest$tree.info <- read.forest$tree.info[subset.idcs,]

  if (!is.null(read.forest$node.obs))
    read.forest$node.obs <- read.forest$node.obs[,subset.idcs]

  return(read.forest)
}

`%<-meta.cache%` <- function(suite, RF.type, verify=c(TRUE, FALSE)) {

  operator <- function(x, new.value) {
    name <- deparse(substitute(x))

    directory <- file.path('assets', suite)
    if (is.null(RF.type) || RF.type == '') {
      filename <- name
    } else {
      filename <- paste(name, RF.type, sep='-')
    }
    base.path <- file.path(directory, filename)
  
    save.path <- paste0(base.path, '.rds')
    if (!file.exists(save.path)) {
      warning(paste(save.path, 'not accessible, regenerating...'))
      if (!dir.exists(directory)) {
        dir.create(directory, recursive=TRUE)
      }
      saveRDS(new.value, save.path)
      assign(name, new.value, inherits=TRUE)
    } else {
      old.value <- readRDS(save.path)
      if (verify) {
        test_that(paste('test if', name, 'is consistent',
                        'for', RF.type), {
          expect_equal(old.value, new.value)
        })
      }
      assign(name, old.value, inherits=TRUE)
    }
  }

  return(operator)
}

make.RF.collection <- function(x, y) {
  `%<-cache%` <- `%<-meta.cache%`('global', NULL, FALSE)

  rand.forest.randomForest %<-cache%
      randomForest::randomForest(Species ~ ., iris)
  
  class.irf <- is.factor(y)
  if (class.irf)
      y <- as.numeric(y) - 1
  rand.forest.ranger %<-cache%
      ranger::ranger(data=cbind(x, y),
                     dependent.variable.name='y',
                     classification=class.irf)
  
  RF.collection <- list(randomForest=rand.forest.randomForest,
                        ranger=rand.forest.ranger)
  return(RF.collection)
}

