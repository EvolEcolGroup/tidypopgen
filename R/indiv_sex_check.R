#' Check the sex of individuals from their sex chromosomes
#'
#' Estimate the sex of each individual from the inbreeding coefficient F of
#' its X chromosome (optionally also using the number of its genotypes on the
#' Y chromosome), and compare it with the sex recorded in the `sex` column of
#' the `gen_tibble`, as PLINK 1.9's `--check-sex`.
#'
#' @details Males have a single X chromosome, so their genotypes on X are
#'   homozygous (stored as homozygotes, as in PLINK), and their F on X is
#'   close to 1; females have two, and an F close to 0. F is computed as in
#'   PLINK 1.9 (`--check-sex`, the same F as `--het`): for each individual,
#'   F = (O - E) / (N - E), where N is the number of loci on X at which it is
#'   genotyped, O the number of those at which it is homozygous, and E the
#'   number expected under Hardy-Weinberg equilibrium (the sum of 1 - 2pq over
#'   those loci, where p is the frequency of the alternate allele, without a
#'   small-sample correction). Loci that are monomorphic (or with no
#'   genotypes) are not used. Allele frequencies are computed as PLINK 1.9
#'   computes them on the X chromosome: from the founders (individuals with
#'   no known parent in the columns `paternal_ID` and `maternal_ID`, if
#'   present), counting one allele for each male (whose heterozygous
#'   genotypes are missing) and two for each other individual.
#'
#'   An individual is male if F is larger than `min_male_xf`, and female if F
#'   is smaller than `max_female_xf`; otherwise (or if it has no genotype on
#'   X), its sex is unknown (NA). With `ycount = TRUE`, the number of its
#'   genotypes (that are not missing) on the Y chromosome must also be at
#'   least `min_male_ycount` for a male, and at most `max_female_ycount` for a
#'   female, as PLINK 1.9's `--check-sex ycount`. Only loci of type "X" and
#'   "Y" are used: loci in the pseudo-autosomal regions should be declared as
#'   such (see [gt_set_chromosome_types()]), as with PLINK's `--split-x`.
#'
#'   In ZW systems (e.g. birds), females are the heterogametic sex: Z and W
#'   take the place of X and Y, and a high F on Z means a female (the
#'   thresholds named after males then apply to females, and vice versa).
#'
#'   The chromosomes must be recognised as sex chromosomes from their names,
#'   or declared with [gt_set_chromosome_types()] (e.g. PLINK's numeric codes,
#'   with `chr_set = "human"`).
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#'   Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC (2007)
#'   PLINK: a tool set for whole-genome association and population-based
#'   linkage analyses. American Journal of Human Genetics, 81(3), 559-575.
#'   \doi{10.1086/519795}
#' @param .x a [`gen_tibble`]; if it has a column `sex` (e.g. as read from a
#'   PLINK .fam file), the sex of each individual is checked against it. The
#'   `sex` column can be a factor or a character vector ("female" and "male",
#'   or "F" and "M"), or PLINK's codes (2 for females, 1 for males, 0 for
#'   unknown).
#' @param max_female_xf the largest F of a female (0.2 by default, as in PLINK
#'   1.9).
#' @param min_male_xf the smallest F of a male (0.8 by default, as in PLINK
#'   1.9).
#' @param ycount boolean, whether to also use the number of genotypes on the Y
#'   chromosome (as PLINK 1.9's `--check-sex ycount`).
#' @param max_female_ycount the largest number of genotypes on Y of a female
#'   (with `ycount = TRUE`; 0 by default, as in PLINK 1.9).
#' @param min_male_ycount the smallest number of genotypes on Y of a male
#'   (with `ycount = TRUE`; 0 by default, as in PLINK 1.9).
#' @param ... currently unused.
#' @returns a tibble with one row per individual (in the order of `.x`) and
#'   the columns of PLINK's `.sexcheck` file: `id`; `sex`, the recorded sex (a
#'   factor with levels "female" and "male"; NA if unknown, or if `.x` has no
#'   `sex` column); `snp_sex`, the sex estimated from the genotypes; `status`,
#'   "OK" if the recorded and estimated sexes are known and the same, and
#'   "PROBLEM" otherwise; `f`, the F of the X chromosome (NaN without
#'   genotypes on X); and, with `ycount = TRUE`, `y_count`, the number of
#'   genotypes on Y that are not missing.
#' @export
#' @examples
#' # three females and two males (hemizygous on X, so homozygous), with
#' # five loci on the X chromosome; the last individual is recorded as a
#' # female, but its genotypes are those of a male
#' x_genotypes <- rbind(
#'   c(1, 0, 2, 1, 1),
#'   c(0, 1, 1, 2, 1),
#'   c(1, 1, 0, 1, 2),
#'   c(2, 0, 2, 0, 2),
#'   c(0, 2, 0, 2, 0)
#' )
#' x_loci <- data.frame(
#'   name = paste0("x", 1:5), chromosome = "X", position = 1:5 * 1000,
#'   genetic_dist = 0, allele_ref = "A", allele_alt = "G"
#' )
#' example_gt <- gen_tibble(x_genotypes,
#'   loci = x_loci,
#'   indiv_meta = data.frame(
#'     id = paste0("ind", 1:5),
#'     sex = c("female", "female", "female", "male", "female")
#'   ),
#'   backingfile = tempfile("sex_check_"),
#'   quiet = TRUE
#' )
#' indiv_sex_check(example_gt)
indiv_sex_check <- function(.x, max_female_xf = 0.2, min_male_xf = 0.8,
                            ycount = FALSE, max_female_ycount = 0,
                            min_male_ycount = 0, ...) {
  rlang::check_dots_empty()
  stopifnot_gen_tibble(.x)
  # the sex of individuals is checked over all of them, whatever their groups
  .x <- dplyr::ungroup(.x)
  # the types of the loci, and the sex chromosomes of the system (XY or ZW)
  locus_types <- loci_type_index(.x$genotypes)
  sex_system <- gt_sex_system(locus_types)
  if (sex_system == "ZW") {
    x_type <- "Z"
    y_type <- "W"
    # the heterogametic sex (1, male in PLINK's codes; 2, female)
    heterogametic <- 2L
  } else {
    x_type <- "X"
    y_type <- "Y"
    heterogametic <- 1L
  }
  x_loci <- which(locus_types == match(x_type, chromosome_type_levels))
  if (length(x_loci) == 0) {
    stop(
      "there are no loci on the X (or Z) chromosome: name it X, or declare ",
      "it with gt_set_chromosome_types()"
    )
  }
  # the recorded sex (1 for males, 2 for females, NA if unknown)
  recorded_sex <- .gt_sex(.x)
  # the founders, whose genotypes give the allele frequencies (as in PLINK,
  # individuals with no known parent)
  founder <- rep(TRUE, nrow(.x))
  for (parent in c("paternal_ID", "maternal_ID")) {
    if (parent %in% names(.x)) {
      founder <- founder & is.na(.x[[parent]])
    }
  }
  # the founders of the heterogametic sex count one allele, the others two
  hemizygous_founder <- founder & (recorded_sex %in% heterogametic)
  diploid_founder <- founder & !(recorded_sex %in% heterogametic)

  fbm <- .gt_get_fbm(.x)
  rows <- .gt_fbm_rows(.x)
  cols <- show_loci(.x)$big_index
  n_indiv <- length(rows)
  # the counts of each individual over the loci on X: genotyped loci (N),
  # heterozygous genotypes, and the expected homozygosity (E)
  n_genotyped <- numeric(n_indiv)
  n_het <- numeric(n_indiv)
  exp_hom <- numeric(n_indiv)
  n_x_used <- 0
  for (block in split_in_blocks(x_loci, bigstatsr::block_size(n_indiv))) {
    # the genotypes of the block (individuals as rows; NA if missing)
    geno <- fbm[rows, cols[block], drop = FALSE]
    # allele frequencies, as PLINK 1.9 on X: hemizygous founders count one
    # allele (a heterozygous genotype is missing), the others two
    hemi <- geno[hemizygous_founder, , drop = FALSE]
    hemi[hemi == 1] <- NA
    dip <- geno[diploid_founder, , drop = FALSE]
    alt_alleles <- colSums(hemi, na.rm = TRUE) / 2 + colSums(dip, na.rm = TRUE)
    n_alleles <- colSums(!is.na(hemi)) + 2 * colSums(!is.na(dip))
    freq <- alt_alleles / n_alleles
    # monomorphic loci (or loci with no genotypes) are not used
    used <- (colSums(!is.na(geno)) > 0) & (n_alleles > 0) &
      (freq >= 1e-8) & (freq <= 1 - 1e-8)
    if (!any(used)) {
      next
    }
    geno <- geno[, used, drop = FALSE]
    freq <- freq[used]
    n_x_used <- n_x_used + sum(used)
    # expected homozygosity of each locus (Nei's, without the small-sample
    # correction, as PLINK 1.9)
    nei <- 1 - 2 * freq * (1 - freq)
    genotyped <- !is.na(geno)
    n_genotyped <- n_genotyped + rowSums(genotyped)
    n_het <- n_het + rowSums(geno == 1, na.rm = TRUE)
    exp_hom <- exp_hom + as.vector(genotyped %*% nei)
  }
  # as PLINK, at least one locus on X must be polymorphic
  if (n_x_used == 0) {
    stop(
      "the sex check requires at least one polymorphic locus on the X (or ",
      "Z) chromosome"
    )
  }
  # F = (O - E) / (N - E), with O = N - heterozygous genotypes
  f <- (n_genotyped - exp_hom - n_het) / (n_genotyped - exp_hom)
  f[n_genotyped == 0] <- NaN
  # the number of genotypes on Y (that are not missing)
  y_count <- NULL
  if (ycount) {
    y_loci <- which(locus_types == match(y_type, chromosome_type_levels))
    y_count <- numeric(n_indiv)
    for (block in split_in_blocks(y_loci, bigstatsr::block_size(n_indiv))) {
      geno <- fbm[rows, cols[block], drop = FALSE]
      y_count <- y_count + rowSums(!is.na(geno))
    }
    if (length(y_loci) == 0) {
      warning("there are no loci on the Y (or W) chromosome")
    }
  }
  # the sex estimated from F (and the counts on Y), with the thresholds of the
  # heterogametic sex (male in PLINK) and of the homogametic sex
  homogametic <- 3L - heterogametic
  snp_sex <- rep(NA_integer_, n_indiv)
  is_hetero <- !is.na(f) & (f > min_male_xf)
  is_homo <- !is.na(f) & !is_hetero & (f < max_female_xf)
  if (ycount) {
    is_hetero <- is_hetero & (y_count >= min_male_ycount)
    is_homo <- is_homo & (y_count <= max_female_ycount)
  }
  snp_sex[is_hetero] <- heterogametic
  snp_sex[is_homo] <- homogametic
  # as PLINK, the status is "OK" only if both sexes are known and the same
  status <- ifelse(
    !is.na(recorded_sex) & !is.na(snp_sex) & (recorded_sex == snp_sex),
    "OK", "PROBLEM"
  )
  sex_factor <- function(codes) {
    factor(c("male", "female")[codes], levels = c("female", "male"))
  }
  res <- tibble::tibble(
    id = .x$id,
    sex = sex_factor(recorded_sex),
    snp_sex = sex_factor(snp_sex),
    status = status,
    f = f
  )
  if (ycount) {
    res$y_count <- as.integer(y_count)
  }
  res
}

# the sex of each individual of a gen_tibble, from its `sex` column: 1 for
# males, 2 for females, and NA if unknown (PLINK's codes, with NA for 0). The
# column can be a factor or a character vector ("male"/"female", "M"/"F", in
# any case, with NA, "", "0", "U" or "unknown" for unknown), or PLINK's codes
# (1, 2, and 0 or NA for unknown); without a sex column, all are unknown
.gt_sex <- function(.x) {
  if (!"sex" %in% names(.x)) {
    return(rep(NA_integer_, nrow(.x)))
  }
  sex <- .x$sex
  if (is.numeric(sex) && !is.factor(sex)) {
    codes <- rep(NA_integer_, length(sex))
    codes[sex %in% 1] <- 1L
    codes[sex %in% 2] <- 2L
    bad <- !is.na(sex) & !(sex %in% c(0, 1, 2))
  } else {
    sex <- tolower(as.character(sex))
    codes <- rep(NA_integer_, length(sex))
    codes[sex %in% c("male", "m", "1")] <- 1L
    codes[sex %in% c("female", "f", "2")] <- 2L
    bad <- !is.na(sex) & !(sex %in% c(
      "male", "m", "1", "female", "f", "2", "", "0", "u", "unknown"
    ))
  }
  if (any(bad)) {
    stop(
      "the 'sex' column should give 'female' or 'male' (or 'F' and 'M', ",
      "or PLINK's codes 2 and 1), not '", .x$sex[bad][1], "'"
    )
  }
  codes
}

# the sex system of the loci of a gen_tibble (from the codes of their types,
# see loci_type_index()): "ZW" if there are loci on Z or W, and "XY"
# otherwise (it is an error to have both)
gt_sex_system <- function(locus_types) {
  present <- chromosome_type_levels[unique(locus_types)]
  xy <- any(c("X", "Y") %in% present)
  zw <- any(c("Z", "W") %in% present)
  if (xy && zw) {
    stop(
      "the gen_tibble has loci on both X/Y and Z/W chromosomes; declare the ",
      "types of chromosomes with gt_set_chromosome_types()"
    )
  }
  if (zw) "ZW" else "XY"
}

# split a vector of indices into blocks of at most block_size elements
split_in_blocks <- function(x, block_size) {
  if (length(x) == 0) {
    return(list())
  }
  split(x, ceiling(seq_along(x) / block_size))
}
