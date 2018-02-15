parallelDuplicateCorrelation.lmer <- function (object, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                                               trim = 0.15, weights = NULL) 
{
  y <- getEAWP(object)
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != narrays) 
    stop("Number of rows of design matrix does not match number of arrays")
  ne <- nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta <- ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    if (any(dim(weights) != dim(M))) 
      weights <- array(weights, dim(M))
    M[weights < 1e-15] <- NA
    weights[weights < 1e-15] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing <- 1
      if (spacing == "rows") 
        spacing <- object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups <- 1
    nspacing <- 1
    Array <- block
  }
  if (is.null(block)) {
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes <- nrow(M)
    if (!is.null(weights)) 
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
  }
  rho <- rep(NA, ngenes)
  nafun <- function(e) NA
  rho = foreach(i = 1:ngenes, .combine = rbind, .packages = c("lme4"), 
                .export = c("M", "Array", "nbeta", "design", "weights", 
                            "nafun"), .verbose = F) %dopar% {
                              y <- drop(M[i, ])
                              o <- is.finite(y)
                              A <- factor(Array[o])
                              nobs <- sum(o)
                              nblocks <- length(levels(A))
                              if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
                                data = cbind(data.frame(y = y[o]),
                                             design[o, , drop = FALSE],
                                             data.frame(A = A))
                                formula1 = paste0('y~',paste(setdiff(colnames(data), c('y','A')), collapse = '+'), '+(1|A)')
                                if (!is.null(weights)) {
                                  w <- drop(weights[i, ])[o]
                                  s <- tryCatch({
                                    mdl = lme4::lmer(formula1, data = data, weights = w);
                                    as.data.frame(VarCorr(mdl))
                                  }, error = nafun)
                                }
                                else {
                                  s <- tryCatch({
                                    mdl = lme4::lmer(formula1, data = data);
                                    as.data.frame(VarCorr(mdl))
                                  }, error = nafun)
                                }
                                if (!is.na(s[1])) 
                                  rho <- s$vcov[1]/sum(s$vcov)
                                else 
                                  rho <- NA
                              }
                            }
  arho <- atanh(pmax(-1, rho))
  mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
  list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
}

parallelDuplicateCorrelation.mm2 <- function (object, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                                              trim = 0.15, weights = NULL) 
{
  y <- getEAWP(object)
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != narrays) 
    stop("Number of rows of design matrix does not match number of arrays")
  ne <- nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta <- ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    if (any(dim(weights) != dim(M))) 
      weights <- array(weights, dim(M))
    M[weights < 1e-15] <- NA
    weights[weights < 1e-15] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing <- 1
      if (spacing == "rows") 
        spacing <- object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups <- 1
    nspacing <- 1
    Array <- block
  }
  if (is.null(block)) {
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes <- nrow(M)
    if (!is.null(weights)) 
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
  }
  rho <- rep(NA, ngenes)
  nafun <- function(e) NA
  rho = foreach(i = 1:ngenes, .combine = rbind, .packages = c("statmod"), 
                .export = c("M", "Array", "nbeta", "design", "weights", 
                            "nafun"), .verbose = F) %dopar% {
                              y <- drop(M[i, ])
                              o <- is.finite(y)
                              A <- factor(Array[o])
                              nobs <- sum(o)
                              nblocks <- length(levels(A))
                              if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 
                                  1) {
                                y <- y[o]
                                X <- design[o, , drop = FALSE]
                                Z <- model.matrix(~0 + A)
                                if (!is.null(weights)) {
                                  w <- drop(weights[i, ])[o]
                                  s <- tryCatch(mixedModel2Fit(y, X, Z, w, only.varcomp = TRUE, 
                                                               maxit = 20)$varcomp, error = nafun)
                                }
                                else {
                                  s <- tryCatch(mixedModel2Fit(y, X, Z, only.varcomp = TRUE, 
                                                               maxit = 20)$varcomp, error = nafun)
                                }
                                if (!is.na(s[1])) 
                                  rho <- s[2]/sum(s)
                                else rho <- NA
                              }
                            }
  arho <- atanh(pmax(-1, rho))
  mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
  list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho)
}

parallelDuplicateCorrelation <- function (object, design = NULL, ndups = 2, spacing = 1, block = NULL, 
                                          trim = 0.15, weights = NULL, method = 'lmer'){
  if (method == 'lmer'){
    parallelDuplicateCorrelation.lmer(object, design, ndups, spacing, block, trim, weights)
  } else if (method == 'mm2'){
    parallelDuplicateCorrelation.mm2(object, design, ndups, spacing, block, trim, weights)
  } else {
    error('Unknown Method')
  }
}