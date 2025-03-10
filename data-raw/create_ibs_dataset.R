# from https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS5.html

offspring.geno <- function(
  n.families,
  n.snps,
  fs = rep(0.5, n.snps),
  n.shared.parents = 2
) {
  #INPUT:
  # n.families, number of families where each family produces two offspring (>0)
  # n.snps, number of independent SNPs used in simulation (>0)
  # fs, vector of allele 1 freqs for SNPs, length == n.snps, values >0 & <1
  # n.shared.parents, 0,1,2 shared parents for the two offspring in each family
  #OUTPUT:
  # X, the genotypes of (2*n.families) offspring, (2*n.families) x n.snps matrix with 0,1,2 entries

  stopifnot(n.families > 0)
  stopifnot(n.snps > 0)
  stopifnot(all(fs > 0 & fs < 1) & length(fs) == n.snps)
  stopifnot(n.shared.parents %in% 0:2)

  if (n.shared.parents == 2) parents = list(c(1, 2), c(1, 2)) #parents[[1]] are the parents of offspring 1
  if (n.shared.parents == 1) parents = list(c(1, 2), c(3, 2))
  if (n.shared.parents == 0) parents = list(c(1, 2), c(3, 4))
  n.parents = 4 - n.shared.parents #4, 3 or 2 for values of n.shared.parent of 0, 1 or 2

  X = matrix(0, nrow = 2 * n.families, ncol = n.snps)
  for (ii in 1:n.families) {
    #each "family" means a pair of offspring that share 'n.shared.parents'
    x.parents = t(replicate(2 * n.parents, rbinom(n.snps, size = 1, prob = fs))) #2*n.parents parental genomes
    for (offs in 1:2) {
      #for two offspring within "family"
      #phase is the indicator of whether offs inherits each parents' 1st allele or not
      phase = t(replicate(2, rbinom(n.snps, size = 1, prob = 0.5))) #phase has one row for each parent
      for (i.parent in 1:2) {
        for (ph in 0:1) {
          loci = (phase[i.parent, ] == ph) #which loci from i.parent have phase ph?
          #add to current offs' genotype i.parent's allele from the correct phase
          X[2 * (ii - 1) + offs, loci] =
            X[2 * (ii - 1) + offs, loci] +
            x.parents[2 * parents[[offs]][i.parent] - ph, loci]
        }
      }
    }
  }
  return(X)
}
set.seed(123)
p = 1000 #SNPs
fs = runif(p, 0.05, 0.5) #MAF at each SNP is Uniform(0.05, 0.5)

X = rbind(
  offspring.geno(n.families = 2, n.snps = p, fs = fs, n.shared.parents = 0),
  offspring.geno(n.families = 2, n.snps = p, fs = fs, n.shared.parents = 1),
  offspring.geno(n.families = 2, n.snps = p, fs = fs, n.shared.parents = 2)
)

X = X[, apply(X, 2, var) > 0] #remove possible monomorphic variants

# sprinkle some missing data
na_n <- 450
na_rows <- sample(dim(X)[1], na_n, replace = TRUE)
na_cols <- sample(dim(X)[2], na_n, replace = TRUE)
for (i in 1:na_n) {
  X[na_rows[i], na_cols[i]] <- NA
}

# transpose it
X = t(X)
fake_bim <- genio::make_bim(n = nrow(X))
fake_fam <- genio::make_fam(n = ncol(X))
# this should be saved directly into inst, I think...
#genio::write_plink("./data-raw/datasets/families", X,bim=fake_bim, fam=fake_fam)
genio::write_plink(
  "./inst/extdata/related/families",
  X,
  bim = fake_bim,
  fam = fake_fam
)
