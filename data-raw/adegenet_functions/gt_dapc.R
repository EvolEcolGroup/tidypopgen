#################
## dapc.genlight
#################
#' @export
gt_dapc <- function(x,
                    pop = NULL,
                    n_pca = NULL,
                    n_da = NULL,
                    scale = FALSE,
                    var_contrib = TRUE,
                    var_loadings = FALSE,
                    pca_info = TRUE,
                    pca_select = c("nbEig", "percVar"),
                    perc_pca = NULL,
                    gt_pca = NULL,
                    ...) {
  pca_select <- match.arg(pca_select)

  if (is.null(pop)) {
    pop.fac <- x$population
  } else {
    pop.fac <- pop
  }

  if (is.null(pop.fac)){
    stop("x does not include pre-defined populations, and `pop' is not provided")
  }

  ## PERFORM PCA ##
  REDUCEDIM <- is.null(gt_pca)

  if (REDUCEDIM) {
    # if no glPca provided
    maxRank <- min(c(nrow(x), length(loci_names(x))))
    pcaX <-
      gt_pca(
        x,
        center = TRUE,
        scale = scale,
        nf = maxRank,
        loadings = FALSE,
        returnDotProd = TRUE,
        ...
      )
  }

  if (!REDUCEDIM) {
    # else use the provided glPca object
    if (is.null(gt_pca$loadings) & var_contrib) {
      warning("Contribution of variables requested but glPca object provided without loadings.")
      var_contrib <- FALSE
    }
    pcaX <- gt_pca
  }

  if (is.null(n_pca)) {
    cumVar <- 100 * cumsum(pcaX$eig) / sum(pcaX$eig)
  }


  ## select the number of retained PC for PCA
  if (!REDUCEDIM) {
    myCol <-
      rep(c("black", "lightgrey"), c(ncol(pcaX$scores), length(pcaX$eig)))
  } else {
    myCol <- "black"
  }

  if (is.null(n_pca) & pca_select == "nbEig") {
    plot(
      cumVar,
      xlab = "Number of retained PCs",
      ylab = "Cumulative variance (%)",
      main = "Variance explained by PCA",
      col = myCol
    )
    cat("Choose the number PCs to retain (>=1): ")
    n_pca <-
      as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  if (is.null(perc_pca) & pca_select == "percVar") {
    plot(
      cumVar,
      xlab = "Number of retained PCs",
      ylab = "Cumulative variance (%)",
      main = "Variance explained by PCA",
      col = myCol
    )
    cat("Choose the percentage of variance to retain (0-100): ")
    nperc.pca <-
      as.numeric(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  ## get n_pca from the % of variance to conserve
  if (!is.null(perc_pca)) {
    n_pca <- min(which(cumVar >= perc_pca))
    if (perc_pca > 99.999)
      n_pca <- length(pcaX$eig)
    if (n_pca < 1)
      n_pca <- 1
  }

  if (!REDUCEDIM) {
    if (n_pca > ncol(pcaX$scores)) {
      n_pca <- ncol(pcaX$scores)
    }
  }


  ## recompute PCA with loadings if needed
  if (REDUCEDIM) {
    pcaX <-
      gt_pca(
        x,
        center = TRUE,
        scale = scale,
        nf = n_pca,
        loadings = var_contrib,
        matDotProd = pcaX$dotProd
      )
  }


  ## keep relevant PCs - stored in XU
  N <- nrow(x)
  X.rank <- sum(pcaX$eig > 1e-14)
  n_pca <- min(X.rank, n_pca)
  if (n_pca >= N)
    n_pca <- N - 1

  U <- pcaX$loadings[, 1:n_pca, drop = FALSE] # principal axes
  XU <- pcaX$scores[, 1:n_pca, drop = FALSE] # principal components
  XU.lambda <-
    sum(pcaX$eig[1:n_pca]) / sum(pcaX$eig) # sum of retained eigenvalues
  names(U) <- paste("PCA-pa", 1:ncol(U), sep = ".")
  names(XU) <- paste("PCA-pc", 1:ncol(XU), sep = ".")


  ## PERFORM DA ##
  ldaX <-
    lda(XU, pop.fac, tol = 1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
  lda.dim <- sum(ldaX$svd ^ 2 > 1e-10)
  ldaX$svd <- ldaX$svd[1:lda.dim]
  ldaX$scaling <- ldaX$scaling[, 1:lda.dim, drop = FALSE]

  if (is.null(n_da)) {
    barplot(
      ldaX$svd ^ 2,
      xlab = "Linear Discriminants",
      ylab = "F-statistic",
      main = "Discriminant analysis eigenvalues",
      col = heat.colors(length(levels(pop.fac)))
    )
    cat("Choose the number discriminant functions to retain (>=1): ")
    n_da <-
      as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  }

  n_da <-
    min(n_da, length(levels(pop.fac)) - 1, n_pca, sum(ldaX$svd > 1e-10)) # can't be more than K-1 disc. func., or more than n.pca
  n_da <- round(n_da)
  predX <- predict(ldaX, dimen = n_da)


  ## BUILD RESULT
  res <- list()
  res$n.pca <- n_pca
  res$n.da <- n_da
  res$tab <- XU
  res$grp <- pop.fac
  res$var <- XU.lambda
  res$eig <- ldaX$svd ^ 2
  res$loadings <- ldaX$scaling[, 1:n_da, drop = FALSE]
  res$means <- ldaX$means
  res$ind.coord <- predX$x
  res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
  res$prior <- ldaX$prior
  res$posterior <- predX$posterior
  res$assign <- predX$class
  res$call <- match.call()


  ## optional: store loadings of variables
  if (pca_info) {
    res$pca.loadings <- as.matrix(U)
    res$pca.cent <- snpbin_list_means(x, alleles_as_units = FALSE)
    if (scale) {
      res$pca.norm <- sqrt(snpbin_list_vars(x, alleles_as_units = FALSE))
    } else {
      res$pca.norm <- rep(1, nLoc(x))
    }
    res$pca.eig <- pcaX$eig
  }

  ## optional: get loadings of variables
  if (var_contrib || var_loadings) {
    var.load <-
      as.matrix(U) %*% as.matrix(ldaX$scaling[, 1:n_da, drop = FALSE])

    if (var_contrib) {
      f1 <- function(x) {
        temp <- sum(x * x)
        if (temp < 1e-12)
          return(rep(0, length(x)))
        return(x * x / temp)
      }
      res$var.contr <- apply(var.load, 2, f1)
    }
    if (var_loadings)
      res$var.load <- var.load
  }

  class(res) <- "dapc"
  return(res)
} # end dapc.genlight
