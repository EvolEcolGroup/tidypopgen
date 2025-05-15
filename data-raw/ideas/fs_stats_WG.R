fs_stats_tg <- function(allele_sharing_mat, pop) {
  Mij <- allele_sharing_mat
  Mii <- diag(Mij) * 2 - 1
  diag(Mij) <- NA
  pop <- factor(pop)
  x <- levels(pop)
  npop <- length(x)
  wil <- lapply(x, function(z) which(pop == z))
  Fi <- lapply(wil, function(x) Mii[x])
  Fsts <- unlist(lapply(wil, function(x) mean(Mij[x, x], na.rm = TRUE)))
  Mb <- 0
  mMij <- matrix(numeric(npop^2), ncol = npop)
  for (i in 2:npop) {
    p1 <- wil[[i]]
    for (j in 1:(i - 1)) {
      p2 <- wil[[j]]
      mMij[i, j] <- mMij[j, i] <- mean(Mij[p1, p2], na.rm = TRUE)
      Mb <- Mb + mMij[i, j]
    }
  }
  diag(mMij) <- Fsts

  # branch off for  pairwise fst
  Fst2x2 <- matrix(NA, ncol = npop, nrow = npop)
  for (i in 2:npop) {
    for (j in 1:(i - 1)) {
      Fst2x2[i, j] <- Fst2x2[j, i] <- ((mMij[i, i] + mMij[j, j]) / 2 - mMij[i, j]) / (1 - mMij[i, j])
    }
  }
  rownames(Fst2x2) <- colnames(Fst2x2) <- levels(pop)


  ######################################
  Mb <- Mb * 2 / (npop * (npop - 1))

  ###################################################



  # branch off for kinship pop pairwise
  resM <- (mMij - Mb) / (1 - Mb)
  rownames(resM) <- colnames(resM) <- levels(pop)

  # branch off for Fis
  for (i in 1:npop) {
    Fi[[i]] <- (Fi[[i]] - Fsts[[i]]) / (1 - Fsts[[i]])
  }
  names(Fi) <- as.character(x)
  Fis <- unlist(lapply(Fi, mean, na.rm = TRUE))
  Fis <- c(Fis, mean(Fis, na.rm = TRUE))

  # branch off for pop fst
  Fsts <- c((Fsts - Mb) / (1 - Mb), mean((Fsts - Mb) / (1 - Mb),
    na.rm = TRUE
  ))
  res <- rbind(Fis, Fsts)
  colnames(res) <- c(levels(pop), "All")
  rownames(res) <- c("Fis", "Fst")
  all.res <- (list(Fi = Fi, FsM = resM, Fst2x2 = Fst2x2, Fs = res))
  class(all.res) <- "fs.dosage"
  all.res
}

fis_stats_tg <- function(allele_sharing_matrix, pop) {
  Mij <- allele_sharing_matrix
  Mii <- diag(Mij) * 2 - 1
  diag(Mij) <- NA
  pop <- factor(pop)
  x <- levels(pop)
  npop <- length(x)
  wil <- lapply(x, function(z) which(pop == z))
  Fi <- lapply(wil, function(x) Mii[x])
  Fsts <- unlist(lapply(wil, function(x) mean(Mij[x, x], na.rm = TRUE)))
  Mb <- 0
  mMij <- matrix(numeric(npop^2), ncol = npop)
  for (i in 2:npop) {
    p1 <- wil[[i]]
    for (j in 1:(i - 1)) {
      p2 <- wil[[j]]
      mMij[i, j] <- mMij[j, i] <- mean(Mij[p1, p2], na.rm = TRUE)
      Mb <- Mb + mMij[i, j]
    }
  }
  diag(mMij) <- Fsts

  Mb <- Mb * 2 / (npop * (npop - 1))

  # estimate individual Fi (inbreeding)
  for (i in 1:npop) {
    Fi[[i]] <- (Fi[[i]] - Fsts[[i]]) / (1 - Fsts[[i]])
  }
  names(Fi) <- as.character(x)
  Fis <- unlist(lapply(Fi, mean, na.rm = TRUE))
  Fis <- c(Fis, mean(Fis, na.rm = TRUE))
  return(Fis)
}



fst_stats_tg <- function(allele_sharing_mat, pop) {
  Mij <- allele_sharing_mat
  Mii <- diag(Mij) * 2 - 1
  diag(Mij) <- NA
  pop <- factor(pop)
  x <- levels(pop)
  npop <- length(x)
  wil <- lapply(x, function(z) which(pop == z))
  Fi <- lapply(wil, function(x) Mii[x])
  Fsts <- unlist(lapply(wil, function(x) mean(Mij[x, x], na.rm = TRUE)))
  Mb <- 0
  mMij <- matrix(numeric(npop^2), ncol = npop)
  for (i in 2:npop) {
    p1 <- wil[[i]]
    for (j in 1:(i - 1)) {
      p2 <- wil[[j]]
      mMij[i, j] <- mMij[j, i] <- mean(Mij[p1, p2], na.rm = TRUE)
      Mb <- Mb + mMij[i, j]
    }
  }
  diag(mMij) <- Fsts

  Mb <- Mb * 2 / (npop * (npop - 1))

  # pop specific fsts
  fst_by_pop <- c((Fsts - Mb) / (1 - Mb), mean((Fsts - Mb) / (1 - Mb),
    na.rm = TRUE
  ))
  return(fst_by_pop)
}
