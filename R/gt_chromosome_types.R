#' Types of chromosomes: autosomes, sex chromosomes and mtDNA
#'
#' The loci of a `gen_tibble` are on autosomes, on sex chromosomes (X and Y
#' in XY systems, Z and W in ZW systems), on their pseudo-autosomal regions
#' (PAR), or on the mitochondrial genome (MT, which can also be used for the
#' genomes of other organelles, e.g. chloroplasts).
#' `gt_set_chromosome_types()` declares the type of chromosomes that can not be
#' recognised from their names (e.g. the numeric codes of PLINK), and the
#' pseudo-autosomal regions of the X chromosome.
#' `gt_chromosome_types()` lists the chromosomes with their type, and
#' `loci_chromosome_types()` gives the type of each locus (e.g. to select the
#' autosomal loci with [select_loci_if()]).
#'
#' @details Chromosomes are recognised from their names, regardless of case
#'   and of a prefix "chr", "chr_", "chromosome" or "chromosome_" (e.g. "X",
#'   "chrX" or "chromosome_x"): "X", "Y", "Z" and "W" are sex chromosomes;
#'   "XY", "PAR", "PAR1" and "PAR2" are pseudo-autosomal regions; and "MT",
#'   "M" and "0M" are the mitochondrial genome. All other chromosomes are
#'   autosomes. Numbers are never recognised, as their meaning depends on the
#'   species (PLINK codes the human X as 23, but chromosome 23 is an autosome
#'   in a dog): declare them with `types` or `chr_set`.
#'
#'   The types declared with `gt_set_chromosome_types()` are kept with the
#'   `gen_tibble` (they are saved by [gt_save()]); each call replaces the
#'   previous declaration, and a call without `types`, `chr_set` and `par`
#'   removes it. A declaration can also bring a recognised chromosome back
#'   among the autosomes (e.g. `types = c(X = "autosome")`).
#'
#'   Functions that pool loci and assume that they are on diploid autosomes
#'   only use the autosomal loci and the pseudo-autosomal regions, as PLINK
#'   does (e.g. for KING), and say how many loci they left out: the
#'   individual statistics [indiv_het_obs()], [indiv_inbreeding()] and
#'   [qc_report_indiv()]; the pairwise statistics [pairwise_king()],
#'   [pairwise_ibs()], [pairwise_dist()], [pairwise_allele_sharing()] and
#'   [pairwise_grm()]; [gt_pca()] (and the deprecated `gt_pca_*SVD()`) and
#'   [gt_pcadapt()]; [gt_extract_f2()], [gt_admixture()] and [gt_snmf()];
#'   [windows_indiv_roh()]; and the population statistics ([pop_fst()],
#'   [pop_fis()], [pop_het_obs()], [pop_het_exp()], [pop_gene_div()],
#'   [pop_global_stats()], [pairwise_pop_fst()], [nwise_pop_pbs()],
#'   [pop_tajimas_d()], the site frequency spectra [pop_sfs()],
#'   [pairwise_pop_sfs()], [nwise_pop_sfs()], [bootstrap_pop_sfs()], and their
#'   `windows_*()` and `jackknife_*()` versions), whose results by locus or by
#'   window then cover the autosomal loci only. [gt_pseudohaploid()] also
#'   tests autosomal loci. The functions that return one value per locus
#'   (`loci_*()`, [qc_report_loci()]), LD clumping and pruning, and
#'   imputation use all loci, without treating the sex chromosomes of males
#'   as haploid. To include a sex chromosome or the mitochondrial genome in
#'   the analyses that use autosomes only, declare it as an autosome.
#'
#'   PLINK codes chromosomes with numbers: with n autosomes, X is n+1, Y is
#'   n+2, the pseudo-autosomal region of X ("XY") is n+3, and MT is n+4 (23
#'   to 26 in humans, as set by `--chr-set n`). `chr_set = n` declares those
#'   codes; `chr_set = "human"` is `chr_set = 22`, plus the codes of
#'   EIGENSOFT for the mitochondrial genome (90) and XY (91).
#'
#'   PLINK's `--split-x` moves the loci of the pseudo-autosomal regions of
#'   the human X chromosome to the chromosome "XY", using the boundaries of a
#'   build of the human genome: a locus on X is in a pseudo-autosomal region
#'   if its position is at most the last position of PAR1, or at least the
#'   first position of PAR2. `par` gives the type "PAR" to those loci without
#'   changing their chromosome: "b36" (or "hg18") uses 2709521 and
#'   154584237, "b37" (or "hg19") 2699520 and 154931044, and "b38" (or
#'   "hg38") 2781479 and 155701383, as PLINK.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#'   Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC (2007)
#'   PLINK: a tool set for whole-genome association and population-based
#'   linkage analyses. American Journal of Human Genetics, 81(3), 559-575.
#'   \doi{10.1086/519795}
#' @param .x a [`gen_tibble`] (for `loci_chromosome_types()`, also its
#'   `genotypes` column).
#' @param types a named character vector giving the type of chromosomes: its
#'   names are the chromosomes (as in [loci_chromosomes()]; for chromosomes
#'   stored as integer codes, their names, see [gt_lookup_chromosomes()]), and
#'   its values their types, one of "autosome", "X", "Y", "Z", "W", "PAR" and
#'   "MT". For example, `c("23" = "X", "24" = "Y")`.
#' @param chr_set the number of autosomes, as in PLINK's `--chr-set`, to
#'   declare the numeric codes of X, Y, XY and MT (see details), or "human".
#'   The types given in `types` take precedence.
#' @param par the pseudo-autosomal regions of the X chromosome: either a build
#'   of the human genome ("b36", "b37", "b38", or their names "hg18", "hg19",
#'   "hg38"), whose boundaries apply to the chromosome of type "X", or a data
#'   frame with columns `chromosome`, `start` and `end` (positions, both
#'   included) giving the regions of each chromosome of type "X" that are
#'   pseudo-autosomal.
#' @returns `gt_set_chromosome_types()` returns the `gen_tibble` with the
#'   declared types; `gt_chromosome_types()` returns a tibble with columns
#'   `chromosome`, `type` and `n_loci` (the number of loci), with one row for
#'   each chromosome and type (an X chromosome with pseudo-autosomal regions
#'   has two rows); `loci_chromosome_types()` returns a factor with the type
#'   of each locus.
#' @rdname gt_chromosome_types
#' @export
#' @examples
#' bed_file <- system.file("extdata", "pop_a.bed", package = "tidypopgen")
#' example_gt <- gen_tibble(bed_file,
#'   backingfile = tempfile("pop_a_"),
#'   quiet = TRUE
#' )
#'
#' # chromosome 23 is taken as an autosome
#' gt_chromosome_types(example_gt)
#'
#' # declare the human codes of PLINK: 23 is the X chromosome
#' example_gt <- gt_set_chromosome_types(example_gt, chr_set = "human")
#' gt_chromosome_types(example_gt)
#' loci_chromosome_types(example_gt)
#'
#' # keep the autosomal loci only
#' example_gt %>%
#'   select_loci_if(loci_chromosome_types(genotypes) == "autosome") %>%
#'   count_loci()
gt_set_chromosome_types <- function(.x, types = NULL, chr_set = NULL,
                                    par = NULL) {
  stopifnot_gen_tibble(.x)
  declared <- character(0)
  # the numeric codes of PLINK (and of EIGENSOFT, for humans)
  if (!is.null(chr_set)) {
    declared <- chr_set_types(chr_set)
  }
  # the types given explicitly take precedence over chr_set
  if (!is.null(types)) {
    types <- check_chromosome_type_names(types)
    declared <- declared[!names(declared) %in% names(types)]
    declared <- c(declared, types)
  }
  # the pseudo-autosomal regions
  if (!is.null(par)) {
    par <- par_regions(par)
  }
  if ((length(declared) == 0) && is.null(par)) {
    attr(.x$genotypes, "chr_types") <- NULL
  } else {
    attr(.x$genotypes, "chr_types") <- list(types = declared, par = par)
  }
  .x
}

#' @rdname gt_chromosome_types
#' @export
gt_chromosome_types <- function(.x) {
  .x <- genotypes_of(.x)
  chromosome <- attr(.x, "loci")$chromosome
  # the chromosome of each locus, as an index into a vector of chromosome
  # names
  chr <- chromosome_index(.x, chromosome)
  locus_types <- loci_type_index(.x)
  # count the loci of each combination of chromosome and type (without a
  # data frame of all loci)
  n_types <- length(chromosome_type_levels)
  combination <- (chr$index - 1L) * n_types + locus_types
  counts <- tabulate(combination, nbins = length(chr$names) * n_types)
  present <- which(counts > 0)
  res <- tibble::tibble(
    chromosome = chr$names[(present - 1L) %/% n_types + 1L],
    type = factor(chromosome_type_levels[(present - 1L) %% n_types + 1L],
      levels = chromosome_type_levels
    ),
    n_loci = counts[present]
  )
  # loci with a missing chromosome
  n_missing <- sum(is.na(chr$index))
  if (n_missing > 0) {
    res <- dplyr::bind_rows(
      res,
      tibble::tibble(
        chromosome = NA_character_,
        type = factor("autosome", levels = chromosome_type_levels),
        n_loci = n_missing
      )
    )
  }
  res
}

#' @rdname gt_chromosome_types
#' @param .col the column to be used when a tibble (or grouped tibble is passed
#'   directly to the function). This defaults to "genotypes" and can only take
#'   that value. There is no need for the user to set it, but it is included
#'   to resolve certain tidyselect operations.
#' @param ... currently unused.
#' @export
loci_chromosome_types <- function(.x, .col = "genotypes", ...) {
  UseMethod("loci_chromosome_types", .x)
}

#' @export
#' @rdname gt_chromosome_types
loci_chromosome_types.tbl_df <- function(.x, .col = "genotypes", ...) {
  stopifnot_has_genotypes(.x)
  .col <- rlang::enquo(.col) %>%
    rlang::quo_get_expr() %>%
    rlang::as_string()
  # confirm that .col is "genotypes"
  if (.col != "genotypes") {
    stop("loci_chromosome_types only works with the genotypes column")
  }
  # extract the column and hand it over to its method
  loci_chromosome_types(.x$genotypes)
}

#' @export
#' @rdname gt_chromosome_types
loci_chromosome_types.vctrs_bigSNP <- function(.x, .col = "genotypes", ...) {
  rlang::check_dots_empty()
  # a factor made from the integer codes of the types (without a vector of
  # strings as long as the loci table)
  structure(loci_type_index(.x),
    levels = chromosome_type_levels,
    class = "factor"
  )
}

# the types of chromosomes, in the order of the levels of the factor returned
# by loci_chromosome_types() (the index of a type in this vector is its code
# in the internal functions below)
chromosome_type_levels <- c("autosome", "X", "Y", "Z", "W", "PAR", "MT")

# the types of the chromosomes that are kept by the functions that only use
# autosomes
autosomal_types <- c("autosome", "PAR")

# the types of chromosomes on which males (or females) are hemizygous, and
# whose haploid genotypes say nothing of the ploidy of an individual
hemizygous_types <- c("X", "Y", "Z", "W", "MT")

# the boundaries of the pseudo-autosomal regions of the human X chromosome in
# each build (the last position of PAR1 and the first position of PAR2), as
# in PLINK's --split-x
par_builds <- list(
  b36 = c(2709521, 154584237),
  b37 = c(2699520, 154931044),
  b38 = c(2781479, 155701383)
)
par_builds$hg18 <- par_builds$b36
par_builds$hg19 <- par_builds$b37
par_builds$hg38 <- par_builds$b38

# check the types given to gt_set_chromosome_types(), and return them as a
# named character vector
check_chromosome_type_names <- function(types) {
  if (!is.character(types) || is.null(names(types)) ||
        anyNA(names(types)) || any(names(types) == "")) {
    stop(
      "'types' should be a named character vector, giving the type ",
      "(the values) of each chromosome (the names)"
    )
  }
  if (anyDuplicated(names(types))) {
    stop("each chromosome should be named only once in 'types'")
  }
  bad <- !types %in% chromosome_type_levels
  if (any(bad)) {
    stop(
      "the types of chromosomes should be one of ",
      paste0("'", chromosome_type_levels, "'", collapse = ", "),
      " (not '", types[bad][1], "')"
    )
  }
  types
}

# the types of the numeric codes of chromosomes set by PLINK's --chr-set
# (with n autosomes, X is n+1, Y n+2, XY n+3 and MT n+4); "human" also
# includes the codes of EIGENSOFT for MT (90) and XY (91)
chr_set_types <- function(chr_set) {
  human <- identical(chr_set, "human")
  if (human) {
    chr_set <- 22
  }
  if (!is.numeric(chr_set) || (length(chr_set) != 1) || is.na(chr_set) ||
        (chr_set < 1) || (chr_set != round(chr_set))) {
    stop(
      "'chr_set' should be the number of autosomes (a positive integer), ",
      "or \"human\""
    )
  }
  types <- stats::setNames(
    c("X", "Y", "PAR", "MT"),
    as.character(chr_set + 1:4)
  )
  if (human) {
    types <- c(types, "90" = "MT", "91" = "PAR")
  }
  types
}

# the pseudo-autosomal regions given to gt_set_chromosome_types(), as a data
# frame with columns chromosome (NA for the chromosome of type X), start and
# end
par_regions <- function(par) {
  if (is.character(par)) {
    if ((length(par) != 1) || !par %in% names(par_builds)) {
      stop(
        "'par' should be one of ",
        paste0("'", names(par_builds), "'", collapse = ", "),
        ", or a data frame with columns 'chromosome', 'start' and 'end'"
      )
    }
    bounds <- par_builds[[par]]
    return(data.frame(
      chromosome = NA_character_,
      start = c(-Inf, bounds[2]),
      end = c(bounds[1], Inf)
    ))
  }
  if (!is.data.frame(par) ||
        !all(c("chromosome", "start", "end") %in% names(par))) {
    stop(
      "'par' should be a build of the human genome, or a data frame with ",
      "columns 'chromosome', 'start' and 'end'"
    )
  }
  if (!is.numeric(par$start) || !is.numeric(par$end) ||
        anyNA(par$start) || anyNA(par$end) || anyNA(par$chromosome) ||
        any(par$start > par$end)) {
    stop(
      "the regions in 'par' should have a chromosome, and numeric 'start' ",
      "and 'end' positions (with 'start' not after 'end')"
    )
  }
  data.frame(
    chromosome = as.character(par$chromosome),
    start = as.numeric(par$start),
    end = as.numeric(par$end)
  )
}

# the codes of the types recognised from the names of chromosomes (and
# declared in chr_types, the declaration of a gen_tibble), one for each name
chromosome_name_type_index <- function(chr_names, chr_types = NULL) {
  types <- rep(1L, length(chr_names))
  # only short names can be those of sex chromosomes or mtDNA (the longest
  # is "chromosome_PAR1"), so that the (slow) regular expression is only
  # applied to a few of the names of, e.g., millions of RAD loci
  short <- which(!is.na(chr_names) & (nchar(chr_names) <= 15))
  if (length(short) > 0) {
    # strip the prefix of the names, as in cast_chromosome_to_int()
    stripped <- toupper(gsub("^(chromosome_|chr_|chromosome|chr)", "",
      chr_names[short],
      ignore.case = TRUE
    ))
    recognised <- c(
      X = "X", Y = "Y", Z = "Z", W = "W",
      XY = "PAR", PAR = "PAR", PAR1 = "PAR", PAR2 = "PAR",
      MT = "MT", M = "MT", "0M" = "MT"
    )
    hit <- stripped %in% names(recognised)
    types[short[hit]] <- match(
      recognised[stripped[hit]],
      chromosome_type_levels
    )
  }
  # the declared types take precedence
  declared <- chr_types$types
  if (length(declared) > 0) {
    hit <- match(chr_names, names(declared))
    types[!is.na(hit)] <- match(
      declared[hit[!is.na(hit)]],
      chromosome_type_levels
    )
  }
  types
}

# the chromosome of each locus, as a list with `index` (for each locus, the
# index of its chromosome in `names`, NA for a missing chromosome) and
# `names` (the names of the chromosomes)
chromosome_index <- function(.x, chromosome) {
  if (is.factor(chromosome)) {
    # a factor already has this structure (a factor used as an index is
    # taken as its integer codes)
    return(list(index = unclass(chromosome), names = levels(chromosome)))
  }
  if (is_int_column(chromosome)) {
    # integer codes: their names are in the chromosome sidecar (if any)
    codes <- unique(chromosome)
    codes <- codes[!is.na(codes)]
    if (file.exists(chromosomes_sidecar_path(.x))) {
      chr_names <- chromosome_sidecar_names(.x, codes)
    } else {
      chr_names <- as.character(codes)
    }
    return(list(index = match(chromosome, codes), names = chr_names))
  }
  chr_names <- unique(as.character(chromosome))
  chr_names <- chr_names[!is.na(chr_names)]
  list(index = match(as.character(chromosome), chr_names), names = chr_names)
}

# the code of the type of each locus of .x (a genotypes column): the type of
# its chromosome, or PAR for the loci of the X chromosome in a
# pseudo-autosomal region
loci_type_index <- function(.x) {
  chr_types <- attr(.x, "chr_types")
  loci <- attr(.x, "loci")
  chromosome <- loci$chromosome
  if (is_int_column(chromosome) && (length(chr_types$types) == 0)) {
    # integer codes without declared types: only the codes of the (few)
    # chromosomes that are not autosomes are looked up (from the chromosome
    # sidecar, read once and cached), as there can be millions of codes
    special <- int_chromosome_types(.x)
    hit <- match(chromosome, special$code)
    types <- rep(1L, length(chromosome))
    types[!is.na(hit)] <- special$type[hit[!is.na(hit)]]
    chr_type_codes <- special$type
    locus_chr_names <- function(i) special$name[hit[i]]
  } else {
    chr <- chromosome_index(.x, chromosome)
    chr_type_codes <- chromosome_name_type_index(chr$names, chr_types)
    types <- chr_type_codes[chr$index]
    # loci with a missing chromosome are taken as autosomal
    types[is.na(types)] <- 1L
    locus_chr_names <- function(i) chr$names[chr$index[i]]
  }
  # loci of the X chromosome in a pseudo-autosomal region
  par <- chr_types$par
  x_code <- match("X", chromosome_type_levels)
  if (!is.null(par) && any(chr_type_codes == x_code)) {
    on_x <- which(types == x_code)
    x_chr <- locus_chr_names(on_x)
    for (i in seq_len(nrow(par))) {
      # a region applies to its chromosome, or to any chromosome of type X
      in_region <- loci$position[on_x] >= par$start[i] &
        loci$position[on_x] <= par$end[i]
      if (!is.na(par$chromosome[i])) {
        in_region <- in_region & (x_chr == par$chromosome[i])
      }
      types[on_x[which(in_region)]] <- match("PAR", chromosome_type_levels)
    }
  }
  types
}

# the integer codes of the chromosomes of .x (a genotypes column with
# integer chromosomes) whose names are recognised as those of chromosomes
# other than autosomes: a data frame with columns code, name and type (the
# code of the type). The names come from the chromosome sidecar, which is
# read once and cached (as long as the file does not change); without a
# sidecar, the names are the codes, and numbers are never recognised
int_chromosome_types <- function(.x) {
  none <- data.frame(code = integer(0), name = character(0), type = integer(0))
  path <- chromosomes_sidecar_path(.x)
  if (!file.exists(path)) {
    return(none)
  }
  info <- file.info(path, extra_cols = FALSE)
  cache_key <- paste0(normalizePath(path), "#types")
  cached <- sidecar_cache[[cache_key]]
  if (!is.null(cached) && identical(cached$size, info$size) &&
        identical(cached$mtime, info$mtime)) {
    return(cached$table)
  }
  sidecar <- data.table::fread(path,
    sep = "\t", header = TRUE, quote = "", na.strings = NULL,
    colClasses = c("integer", "character"), showProgress = FALSE,
    data.table = FALSE
  )
  types <- chromosome_name_type_index(sidecar[[2]])
  special <- types != 1L
  res <- data.frame(
    code = sidecar[[1]][special],
    name = sidecar[[2]][special],
    type = types[special]
  )
  assign(cache_key,
    list(size = info$size, mtime = info$mtime, table = res),
    envir = sidecar_cache
  )
  res
}

# for the loci of .x (a genotypes column), whether each is autosomal (or in
# a pseudo-autosomal region), or NULL if they all are. Most gen_tibbles have
# no chromosome of another type, which is checked first on the levels of the
# chromosomes (or on the cached codes of integer chromosomes), so that the
# type of each locus is only found when it is needed
autosomal_loci <- function(.x) {
  chr_types <- attr(.x, "chr_types")
  chromosome <- attr(.x, "loci")$chromosome
  keep_codes <- match(autosomal_types, chromosome_type_levels)
  if (is.factor(chromosome)) {
    level_types <- chromosome_name_type_index(levels(chromosome), chr_types)
    other <- !level_types %in% keep_codes
    # (levels of other types may be unused, e.g. after selecting loci)
    if (!any(other) ||
          !any(tabulate(chromosome, nbins = nlevels(chromosome))[other] > 0)) {
      return(NULL)
    }
  } else if (is_int_column(chromosome) && (length(chr_types$types) == 0)) {
    special <- int_chromosome_types(.x)
    if (!any(chromosome %in% special$code[!special$type %in% keep_codes])) {
      return(NULL)
    }
  }
  keep <- loci_type_index(.x) %in% keep_codes
  if (all(keep)) {
    return(NULL)
  }
  keep
}

# the gen_tibble (or genotypes column) .x restricted to its autosomal loci
# (and those in pseudo-autosomal regions), with a message giving the number
# of loci that are left out, for the functions that assume that loci are on
# diploid autosomes (as PLINK does for, e.g., KING); .x is returned as it is
# if all its loci are autosomal (or if it is not a gen_tibble, e.g. the
# name of a file, which is left to the caller)
.gt_autosomes_only <- function(.x) {
  if (inherits(.x, "vctrs_bigSNP")) {
    geno <- .x
  } else if (is.data.frame(.x) && inherits(.x$genotypes, "vctrs_bigSNP")) {
    geno <- .x$genotypes
  } else {
    return(.x)
  }
  keep <- autosomal_loci(geno)
  if (is.null(keep)) {
    return(.x)
  }
  if (!any(keep)) {
    stop(
      "there are no autosomal loci (this function only uses autosomes): ",
      "declare the chromosomes to use as autosomes with ",
      "gt_set_chromosome_types()"
    )
  }
  left_out <- loci_type_index(geno)[!keep]
  message(
    "Using autosomal loci only: ", length(left_out), " loci on ",
    paste(chromosome_type_levels[sort(unique(left_out))], collapse = ", "),
    " are left out (see ?gt_chromosome_types)"
  )
  loci <- vctrs::vec_slice(attr(geno, "loci"), which(keep))
  if (inherits(.x, "vctrs_bigSNP")) {
    attr(.x, "loci") <- loci
  } else {
    attr(.x$genotypes, "loci") <- loci
  }
  .x
}

# the declared types of chromosomes of two gen_tibbles that are bound
# together (by rbind()): their union, which must not give two types to a
# chromosome, nor two different sets of pseudo-autosomal regions
merge_chr_types <- function(chr_types_1, chr_types_2) {
  if (is.null(chr_types_1)) {
    return(chr_types_2)
  }
  if (is.null(chr_types_2)) {
    return(chr_types_1)
  }
  types_1 <- chr_types_1$types
  types_2 <- chr_types_2$types
  shared <- intersect(names(types_1), names(types_2))
  if (any(types_1[shared] != types_2[shared])) {
    stop(
      "the gen_tibbles declare different types for some chromosomes ",
      "(see gt_set_chromosome_types())"
    )
  }
  par_1 <- chr_types_1$par
  par_2 <- chr_types_2$par
  if (!is.null(par_1) && !is.null(par_2) && !identical(par_1, par_2)) {
    stop(
      "the gen_tibbles declare different pseudo-autosomal regions ",
      "(see gt_set_chromosome_types())"
    )
  }
  list(
    types = c(types_1, types_2[setdiff(names(types_2), shared)]),
    par = if (is.null(par_1)) par_2 else par_1
  )
}
