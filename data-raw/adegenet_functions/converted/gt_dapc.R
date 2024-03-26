#' Discriminant Analysis of Principal Components for gen_tibble
#'
#' This function implements the Discriminant Analysis of Principal Components
#' (DAPC, Jombart et al. 2010). This method descibes the diversity between
#' pre-defined groups. When groups are unknown, use [gt_pca_find_clusters()] to
#' infer genetic clusters. See 'details' section for a succint
#' description of the method, and the vignette in the package `adegenet`
#' ("adegenet-dapc") for a
#' tutorial. This function returns objects of class [`adegenet::dapc`] which are compatible
#' with methods from `adegenet`; graphical methods for DAPC are documented in
#' [adegenet::scatter.dapc]
#' (see ?scatter.dapc).
#'
#' The Discriminant Analysis of Principal Components (DAPC) is designed to
#' investigate the genetic structure of biological populations. This
#' multivariate method consists in a two-steps procedure. First, genetic
#' data are transformed (centred, possibly scaled) and submitted to a
#' Principal Component Analysis (PCA). Second, principal components of
#' PCA are submitted to a Linear Discriminant Analysis (LDA). A trivial matrix
#' operation allows to express discriminant functions as linear combination
#' of alleles, therefore allowing one to compute allele contributions. More
#' details about the computation of DAPC are to be found in the indicated reference.
#'
#' @references Jombart T, Devillard S and Balloux F (2010) Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. BMC Genetics 11:94. doi:10.1186/1471-2156-11-94
#'
#'
#' @param x an object of class `gt_pca`, or its subclass `gt_pca_clust`
#' @param pop either a factor indicating the group membership of individuals;
#' or an integer defining the desired *k* if x is a `gt_pca_clust`; or NULL, if
#' 'x' is a `gt_pca_clust` and contain an element 'best_k',
#' usually generated with [gt_pca_best_k()],
#' which will be used to select the clustering level.
#' @param n_pca number of principal components to be used in the Discriminant
#' Analysis. If NULL, k-1 will be used.
#' @param n_da an integer indicating the number of axes retained in the
#' Discriminant Analysis step.
#' @param var_contrib a logical indicating whether the contribution of
#' loci should be stored
#' (TRUE, default) or not (FALSE). Such output can be useful, but can
#' also create huge matrices when there are a lot of loci.
#' @param var_loadings a logical indicating whether the loadings of loci should
#' be stored (TRUE) or
#' not (FALSE, default). Such output can be useful, but can also create
#' huge matrices when there are a lot of loci.
#' @param pca_info a logical indicating whether information about the prior
#' PCA should be stored (TRUE, default) or not (FALSE). This information
#' is required to predict group membership of new individuals using predict,
#' but makes the object slightly bigger.
#' @returns an object of class [adegenet::dapc]
#' @export

# AM: Thoughts about data structures. The original DAPC blended pca info
# within the object. For a cleaner job at predicting, it would be best to
# simply store the pca object as an element within the object. This would break
# compatibility, but only with the predict function, which in any case will not
# work.
# Once we have tidiers and autoplot, it might be best to reorganise the
# object to fully fit our purposes, and give up trying to be backcompatible
# with adegenet


gt_dapc <- function(x, pop = NULL, n_pca = NULL, n_da=NULL,
                          var_contrib=TRUE,
                      var_loadings=FALSE, pca_info =TRUE){
  if (!inherits(x,"gt_pca")){
    stop("'x' should be a 'gt_pca' object")
  }
  if (is.null(x$center)){
    stop("'x' was run without centering; centering is necessary for 'gt_dapc'")
  }

  if(is.null(pop)) { # if no pop was given, use best_k
    if (any(!inherits(x,"gt_pca_clust"),is.null(x$best_k))){
      stop("if 'pop' is not set, 'x' should be a 'gt_pca_clust ")
    }
    pop.fac <- as.factor(x$clusters$groups[[x$best_k]])
  } else if (is.factor(pop)) { # if a factor with all assignments was given
    pop.fac <- pop
  } else if (is.numeric(pop) & inherits(x,"gt_pca_clust")) { # if pop is the k value
    pop.fac <- as.factor(x$clusters$groups[[pop]])
  }

  if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")
  n_pop <- nlevels(pop)


  if (is.null(n_pca) ){
    if (inherits(x,"gt_pca_clust")){ # if we generated clusters, use the same pca
      n_pca   <- x$clusters$n_pca
    } else { # use all principal components
      n_pca   <- length(x$d)
    }
    # by default, use number of clusters minus 1
    if (n_pca>n_pop){
      n_pca <- n_pop - 1
    }
  } else { # if n_pca was given, check that it is not too large
    if(n_pca > ncol(x$u)) {
      n_pca <- ncol(x$u)
    }
  }


    if(is.null(x$v) & var_contrib) {
      warning("'var_contrib' is set to true, but 'x' object provided without loadings ($v).")
      var_contrib <- FALSE
    }


  U <- x$u[, 1:n_pca ,  drop=FALSE] # principal axes
  XU <- sweep(x$u, 2, x$d, '*')[, 1:n_pca ,  drop=FALSE] # principal components
  # note taht this is the proportion of variance out of the variance we started with (i.e. what we retained with the PCAs)
  XU.lambda <- sum(x$d[1:n_pca ] )/sum(x$d) # sum of retained eigenvalues
  names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
  names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


  ## PERFORM DA ##
  ldaX <- MASS::lda(XU, pop.fac, tol=1e-30) # tol=1e-30 is a kludge, but a safe (?) one to avoid fancy rescaling by lda.default
  lda.dim <- sum(ldaX$svd^2 > 1e-10)
  ldaX$svd <- ldaX$svd[1:lda.dim]
  ldaX$scaling <- ldaX$scaling[,1:lda.dim,drop=FALSE]

  # @TODO we need to bring this back as an autoplot for this object type.
  # if(is.null(n.da)){
  #   barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(pop.fac))) )
  #   cat("Choose the number discriminant functions to retain (>=1): ")
  #   n.da <- as.integer(readLines(con = getOption('adegenet.testcon'), n = 1))
  # }

  n_da <- min(n_da, length(levels(pop.fac ))-1, n_pca ,  sum(ldaX$svd>1e-10)) # can't be more than K-1 disc. func., or more than n.pca
  n_da <- round(n_da)
  predX <- stats::predict(ldaX, dimen=n_da)


  ## BUILD RESULT
  res <- list()
  res$n.pca <- n_pca
  res$n.da <- n_da
  res$tab <- XU
  res$grp <- pop.fac
  res$var <- XU.lambda # note that this is the variance out fo the variance that was captured by the retained PCs
  res$eig <- ldaX$svd^2
  res$loadings <- ldaX$scaling[, 1:n_da, drop=FALSE]
  res$means <- ldaX$means
  res$ind.coord <-predX$x
  res$grp.coord <- apply(res$ind.coord, 2, tapply, pop.fac, mean)
  res$prior <- ldaX$prior
  res$posterior <- predX$posterior
  res$assign <- predX$class
  res$call <- match.call()

  # ## optional: store loadings of variables
  # @BUG we need to sort out the slots as these are not correct
  # @TODO our objects are missing several of these slots
   if(pca_info){
     warning("conversion of objects slots is inconmplete, don't use this option yet!")
    res$pca.loadings <- as.matrix(U)
     # res$pca.cent <- x$cent
     # if(!is.null(x$norm)) {
     #   res$pca.norm <- x$norm
     # } else {
     #   res$pca.norm <- rep(1, length(x$cent))
     # }
     res$pca.eig <- x$d^2 # TODO check, this should get back the eigen from glPCA
     # note that the default allele.as.unit is FALSE for glPCA
  }

  ## optional: get loadings of variables
  if(var_contrib || var_loadings){
    message("conversion of objects slots needs to be tested for this option")
    var.load <- as.matrix(U) %*% as.matrix(ldaX$scaling[,1:n_da,drop=FALSE])

    if(var_contrib){
      f1 <- function(x){
        temp <- sum(x*x)
        if(temp < 1e-12) return(rep(0, length(x)))
        return(x*x / temp)
      }
      res$var.contr <- apply(var.load, 2, f1)
    }
    if(var_loadings) res$var.load <- var.load
  }

  class(res) <- "dapc"
  return(res)
}

