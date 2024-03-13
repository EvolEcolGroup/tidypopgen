raw_path_pop_a <- system.file("extdata/pop_a.bed", package = "tidypopgen")
bigsnp_path_a <- bigsnpr::snp_readBed(raw_path_pop_a, backingfile = tempfile("test_a_"))
pop_a_snp <- bigsnpr::snp_attach(bigsnp_path_a)

geno_start <- pop_a_snp$genotypes[]

# quick function to swap loci
swap_locus <- function(x){
  x_new <- x
  x_new[x==2]<-0
  x_new[x==0]<-2
  x_new
}

test_that("subset_bigSNP correctly subsets the data", {
  indiv_indices=c(2,1,4)
  loci_indices=c(1:5,6,11,10)
  swap_indices = c(2,5)
  new_snp <- subset_bigSNP(pop_a_snp, indiv_indices = indiv_indices,
                           loci_indices = loci_indices,
                           swap_indices = swap_indices)
  # we first swap the old loci
  geno_start[,swap_indices] <- swap_locus(geno_start[,swap_indices])
  expect_identical(new_snp$genotypes[],geno_start[indiv_indices,loci_indices])
  #check the swapped locus map
  swapped_locus1 <- pop_a_snp$map$marker.ID[swap_indices[1]]
  expect_true(pop_a_snp$map %>% dplyr::filter(marker.ID==swapped_locus1) %>% pull(allele1)==
    new_snp$map %>% dplyr::filter(marker.ID==swapped_locus1) %>% pull(allele2))

})




