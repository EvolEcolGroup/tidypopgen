# Time each step of gen_tibble() on a PLINK .bed file (the steps of
# gen_tibble.character() and gen_tibble_bed_rds(), with timers)
# usage: Rscript breakdown_full.R <lib> <bed file> <out dir> <n_cores> <names_as_int>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(tidypopgen, lib.loc = args[1]))
ns <- asNamespace("tidypopgen")
bed <- args[2]
out <- file.path(args[3], "full")
n_cores <- as.integer(args[4])
names_as_int <- as.logical(args[5])
bim <- sub("\\.bed$", ".bim", bed)
fam <- sub("\\.bed$", ".fam", bed)
tm <- function(label, expr) {
  t0 <- proc.time()
  res <- force(expr)
  t1 <- proc.time() - t0
  cat(sprintf("%-30s elapsed %7.2f  user %7.2f  sys %7.2f  (%s)\n",
              label, t1[["elapsed"]], t1[["user.self"]], t1[["sys.self"]],
              format(Sys.time(), "%T")))
  invisible(res)
}
tmp <- bigstatsr::FBM.code256(2, 2, code = bigsnpr::CODE_012,
                              backingfile = tempfile(tmpdir = args[3]))
unlink(paste0(out, c(".bk", ".rds", ".gt")))
t_all <- proc.time()
if (names_as_int) {
  bim_df <- tm("fread .bim", ns$read_loci_file(bim, character_cols = c(5, 6),
                                                drop = 2))
  names(bim_df)[1:5] <- c("chromosome", "genetic.dist", "physical.pos",
                          "allele1", "allele2")
} else {
  bim_df <- tm("fread .bim", ns$read_loci_file(bim,
                                                character_cols = c(2, 5, 6)))
  names(bim_df)[1:6] <- c("chromosome", "marker.ID", "genetic.dist",
                          "physical.pos", "allele1", "allele2")
}
loci <- tm("loci_from_bim (rest)",
           ns$loci_from_bim(bim_df, names_as_int = names_as_int))
rm(bim_df)
meta <- tm("indiv_meta_from_fam", ns$indiv_meta_from_fam(fam))
has_n_cores <- "n_cores" %in% names(formals(ns$fbm_read_bed))
fbm_path <- tm("fbm_read_bed", if (has_n_cores) {
  ns$fbm_read_bed(bed, n_indiv = nrow(meta), n_snp = nrow(loci),
                  backingfile = out, n_cores = n_cores)
} else {
  ns$fbm_read_bed(bed, n_indiv = nrow(meta), n_snp = nrow(loci),
                  backingfile = out)
})
x <- tm("build gen_tibble", {
  fbm_obj <- readRDS(fbm_path)
  meta$genotypes <- ns$new_vctrs_bigsnp(
    fbm_obj = fbm_obj, fbm_file = fbm_path, loci = loci,
    indiv_id = meta$id, ploidy = 2, fbm_ploidy = rep(2, nrow(meta))
  )
  tibble::new_tibble(meta, class = "gen_tbl")
})
tm("remove_stale_sidecars", ns$remove_stale_sidecars(gt_get_file_names(x)[2]))
tm("has_genotypes_without_alleles", ns$has_genotypes_without_alleles(x))
tm("gt_save_light", ns$gt_save_light(x, quiet = TRUE))
t_all <- proc.time() - t_all
cat(sprintf("TOTAL elapsed %.2f user %.2f sys %.2f; loci %d; R memory %.1f GB\n",
            t_all[["elapsed"]], t_all[["user.self"]], t_all[["sys.self"]],
            nrow(loci), sum(gc()[, 2]) / 1024))
rm(x, fbm_obj)
unlink(paste0(out, c(".bk", ".rds", ".gt")))
unlink(tmp$backingfile)
