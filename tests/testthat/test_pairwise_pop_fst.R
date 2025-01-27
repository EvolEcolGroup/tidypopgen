library(KRIS)

test_genotypes <- rbind(c(1,1,0,1,1,0),
                        c(2,1,0,NA,0,0),
                        c(2,NA,0,0,1,1),
                        c(1,0,0,1,0,0),
                        c(1,2,0,1,2,1),
                        c(0,0,0,0,NA,1),
                        c(0,1,1,0,1,NA))
test_indiv_meta <- data.frame (id=c("a","b","c","d","e","f","g"),
                               population = c("pop1","pop1","pop2","pop2","pop1","pop3","pop3"))
test_loci <- data.frame(name=paste0("rs",1:6),
                        chromosome=paste0("chr",c(1,1,1,1,2,2)),
                        position=as.integer(c(3,5,65,343,23,456)),
                        genetic_dist = as.double(rep(0,6)),
                        allele_ref = c("A","T","C","G","C","T"),
                        allele_alt = c("T","C", NA,"C","G","A"))

test_gt <- gen_tibble(x = test_genotypes, loci = test_loci, indiv_meta = test_indiv_meta, quiet = TRUE)


test_that("pairwise_pop_fst Nei 87 computes correctly",{
  test_gt <- test_gt %>% dplyr::group_by(population)
  test_hier <- gt_as_hierfstat(test_gt)
  # compare results against hierfstat for Nei87 (Nei86 does not correct for Ho
  # when computing Ht, so it gives a different result)
  nei_gt <- test_gt %>% pairwise_pop_fst(method="Nei87")

  nei_hier <- hierfstat::pairwise.neifst(test_hier)
  # hiefstat values are rounded to 4 dp
  expect_true(all.equal(tidy_dist_matrix(nei_hier)$value, round(nei_gt$value,4)))

})

test_that("pairwise_pop_fst Hudson method computes correctly",{
  # We use the fst.hudson function from KRIS to compute pairwise Fst
  # This is a version with the ratio of means; it is meant to be the default version
  # in KRIS, but there is a problem with a double declaration of the function and the
  # wrong version is used for the package in CRAN. KRIS is now archived, so the CRAN
  # package can not be fixed

  fst.hudson <-function(X, idx.p1, idx.p2){
    prestep.fst.one.marker <- function(alleles,idx.p1,idx.p2){

      #Pop 1
      G = alleles[idx.p1]
      no.AA=length(which(G==0))
      no.AB=length(which(G==1))
      no.BB=length(which(G==2))
      n1=(no.AA+no.AB+no.BB)*2
      p.A=(no.AA*2 + no.AB)/n1
      #p.B=(no.BB*2 + no.AB)/n1
      #p1 = min(p.A,p.B,na.rm = T)
      p1 = p.A

      #Pop 2
      G = alleles[idx.p2]
      no.AA=length(which(G==0))
      no.AB=length(which(G==1))
      no.BB=length(which(G==2))
      n2=(no.AA+no.AB+no.BB)*2
      p.A=(no.AA*2 + no.AB)/n2
      #p.B=(no.BB*2 + no.AB)/n2
      #p2 = min(p.A,p.B,na.rm = T)
      p2 = p.A


      N = (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
      D = p1*(1-p2) + p2*(1-p1)

      return(c(N,D))
    }

    set.fst = apply(X,2,prestep.fst.one.marker,idx.p1=idx.p1, idx.p2=idx.p2)

    # Bhatia, et al 2013, "This is the basis of our recommendation that FST be estimated as a ratio of averages."

    fst = mean(set.fst[1,],na.rm=T) / mean(set.fst[2,],na.rm=T)
    #fst = mean(set.fst[1,]/set.fst[2,],na.rm = T)

    return(fst)
  }


  # Load anolis data: 46 individuals and 3249 variants
  vcf_path <- system.file("/extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
                          package = "tidypopgen")
  anole_gt <- gen_tibble(vcf_path, quiet = TRUE, backingfile = tempfile("anolis_"))
  # attach metadata
  pops_path <- system.file("/extdata/anolis/plot_order_punctatus_n46.csv",
                           package = "tidypopgen")
  pops <- readr::read_csv(pops_path)
  anole_gt <- anole_gt %>% mutate(id = gsub('punc_',"",.data$id,))
  anole_gt <- anole_gt %>% mutate(population = pops$pop[match(pops$ID,.data$id)])
  # group by population
  anole_gt <- anole_gt %>% dplyr::group_by(population)
  # impute
  anole_gt <- anole_gt %>% gt_impute_simple(method = "mode")
  gt_set_imputed(anole_gt, TRUE)
  # calculate pairwise Fst in tidypopgen
  hudson <- anole_gt %>% pairwise_pop_fst(method="Hudson")

  # select individuals from each population
  af <- which(anole_gt$population == 'AF')
  eam <- which(anole_gt$population == 'Eam')
  wam <- which(anole_gt$population == 'Wam')
  # get genotypes matrix
  genotypes <- show_genotypes(anole_gt)
  # calculate pairwise Fst in KRIS
  kris_af_eam <- fst.hudson(genotypes,af,eam)
  kris_af_wam <- fst.hudson(genotypes,af,wam)
  kris_eam_wam <- fst.hudson(genotypes,eam,wam)


  expect_equal(hudson$value[1],kris_af_eam)
  expect_equal(hudson$value[2],kris_af_wam)
  expect_equal(hudson$value[3],kris_eam_wam)

  # And by locus
  # tidypop_locus_fst <- anole_gt %>% pairwise_pop_fst(method = "Hudson",by_locus = TRUE)
  # kris_locus_fst <- fst.each.snp.hudson(genotypes,af,eam)
  # kris_locus_fst <- fst.each.snp.hudson(genotypes,af,wam)
  # kris_locus_fst <- fst.each.snp.hudson(genotypes,eam,wam)

})


