#' Look up the original names of loci and chromosomes stored as integers
#'
#' With `names_as_int = TRUE` and `chromosomes_as_int = TRUE`, [gen_tibble()]
#' stores the names of the loci, and their chromosomes, as integers, to save
#' memory and time with large datasets (e.g. RAD loci assembled de novo, with
#' millions of loci, each its own "chromosome"). The original names are
#' written to sidecar files next to the backing file (see Details), from which
#' these functions get them back:
#' - `gt_lookup_loci()` gives back the original names and chromosomes of some
#'   loci (e.g. the hits of a selection scan), without changing the
#'   `gen_tibble`;
#' - `gt_lookup_chromosomes()` gives back the names of chromosome codes (e.g.
#'   the `chromosome` column of the results of window statistics or runs of
#'   homozygosity);
#' - `gt_restore_loci()` puts the original names back into the loci table of
#'   a `gen_tibble` (e.g. once it has been filtered down to fewer loci).
#'
#' @details The sidecars of the backing file `foo.bk` are:
#' - `foo.loci.tsv` (written with `names_as_int = TRUE`), with the columns
#'   `int_name`, `name`, `chromosome` and `position`, and one line per locus,
#'   in increasing order of integer name. The names are those that
#'   `names_as_int = FALSE` would give (with a missing ID replaced by the
#'   chromosome and position), and the chromosomes are their names (not their
#'   codes);
#' - `foo.chromosomes.tsv` (written with `chromosomes_as_int = TRUE`), with
#'   the columns `int_chromosome` and `chromosome`, and one line per
#'   chromosome, in increasing order of code.
#'
#' Both are tab separated text files with a header, with no quotes, and "."
#' for missing values, so they can be read with any tool (e.g.
#' `data.table::fread(path, quote = "")` or `readr::read_tsv(path, quote =
#' "")`), and joined to [show_loci()] by `name = int_name` or `chromosome =
#' int_chromosome`.
#'
#' The sidecars are written by [gen_tibble()] when it reads a VCF file with
#' `parser = "cpp"`, and by [gt_merge_vcfs()] and [gt_concat_vcfs()].
#' Functions that write a new backing file ([gt_update_backingfile()],
#' [gt_order_loci()], [rbind()][rbind.gen_tbl]) write new sidecars for it,
#' with the loci and chromosomes that they keep (integer names and codes keep
#' their values, so they still find their line). Keep the sidecars with the
#' other files of the `gen_tibble` (see [gt_save()]): without them, the
#' original names are lost.
#'
#' Only the file name ties a sidecar to its data, so the position of each
#' locus in the locus sidecar is compared with its position in the loci
#' table: a different position means that the sidecar belongs to other data
#' (e.g. files copied by hand from different runs), and is an error. If you
#' changed the positions of the loci yourself (e.g. after a liftover to a new
#' assembly), the names in the sidecar are still right: use `check_positions
#' = FALSE`. The chromosome sidecar has no positions to check.
#'
#' Lookups are fast: the first lookup in a session reads through a sidecar
#' once to index it (about a second per GB), and later lookups of a few loci
#' only read the parts of the file that hold them.
#'
#' @param .x a [`gen_tibble`] (or, except for `gt_restore_loci()`, its
#'   `genotypes` column)
#' @param loci a subset of the loci table of `.x`: rows of [show_loci()], in
#'   any order, possibly with extra columns (e.g. the statistics of a scan);
#'   or row numbers of `show_loci(.x)`, as given by [which()].
#' @param chromosomes for `gt_lookup_chromosomes()`, a vector of chromosome
#'   codes (e.g. the `chromosome` column of the results of
#'   [windows_pairwise_pop_fst()] or [gt_roh_window()]); for
#'   `gt_restore_loci()`, boolean: whether to restore the names of the
#'   chromosomes (if they are integer codes), default TRUE.
#' @param names boolean: whether to restore the names of the loci (if they
#'   are integers). Default is TRUE.
#' @param check_positions boolean: whether to check that the positions of the
#'   loci in the locus sidecar are the same as in `loci` (or in the loci table
#'   of `.x`, for `gt_restore_loci()`). For `gt_lookup_loci()`, they are only
#'   checked if `loci` has a `position` column. Default is TRUE.
#' @returns `gt_lookup_loci()` returns `loci` (as a table), with integer names
#'   in its `name` column, and integer codes in its `chromosome` column,
#'   replaced by their original names (as character vectors); the other
#'   columns, and the rows and their order, are unchanged. Columns that
#'   already hold names are left as they are, so the same code works on data
#'   read with or without the integer options.
#'
#'   `gt_lookup_chromosomes()` returns the names of the chromosomes (a
#'   character vector, with NA for NA). Chromosomes that are already names (a
#'   factor or a character vector) are returned as a character vector.
#'
#'   `gt_restore_loci()` returns `.x`, with the original names in its loci
#'   table (as if it had been read without the integer options: names as a
#'   character vector, and chromosomes as a factor). It does not save it: use
#'   [gt_save()].
#' @rdname gt_lookup_loci
#' @export
#' @examples
#' vcf_path <- system.file(
#'   "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
#'   package = "tidypopgen"
#' )
#' anole_gt <- gen_tibble(vcf_path,
#'   backingfile = tempfile("anolis_"),
#'   names_as_int = TRUE, chromosomes_as_int = TRUE, quiet = TRUE
#' )
#' head(show_loci(anole_gt))
#'
#' # the original names of the loci with the most missing genotypes
#' missing <- loci_missingness(anole_gt)
#' gt_lookup_loci(anole_gt, which(missing > 0.6))
#'
#' # the same, from a subset of the loci table, with the missingness
#' hits <- show_loci(anole_gt) %>%
#'   dplyr::mutate(missing = missing) %>%
#'   dplyr::filter(missing > 0.6)
#' gt_lookup_loci(anole_gt, hits)
#'
#' # the names of some chromosome codes
#' gt_lookup_chromosomes(anole_gt, c(1, 5, 5, NA))
#'
#' # put the original names back into the gen_tibble
#' anole_gt <- gt_restore_loci(anole_gt)
#' head(show_loci(anole_gt))
gt_lookup_loci <- function(.x, loci, check_positions = TRUE) {
  .x <- genotypes_of(.x)
  if (is.data.frame(loci)) {
    # used as it is
  } else if (is.numeric(loci)) {
    # row numbers of the loci table
    n_loci <- nrow(show_loci(.x))
    if (anyNA(loci) || any(loci < 1 | loci > n_loci) ||
          any(loci != round(loci))) {
      stop(
        "'loci' should be row numbers of show_loci(.x), between 1 and ",
        n_loci
      )
    }
    loci <- show_loci(.x)[loci, ]
  } else {
    stop(
      "'loci' should be a subset of the loci table of .x (a data frame), ",
      "or row numbers of show_loci(.x)"
    )
  }
  if (!any(c("name", "chromosome") %in% names(loci))) {
    stop(
      "'loci' should have a 'name' or 'chromosome' column, as the loci ",
      "table of .x"
    )
  }
  if (is_int_column(loci$name)) {
    positions <- NULL
    if (check_positions && ("position" %in% names(loci))) {
      positions <- loci$position
    }
    loci$name <- loci_sidecar_names(.x, loci$name, positions)
  }
  if (is_int_column(loci$chromosome)) {
    loci$chromosome <- chromosome_sidecar_names(.x, loci$chromosome)
  }
  loci
}

#' @rdname gt_lookup_loci
#' @export
gt_lookup_chromosomes <- function(.x, chromosomes) {
  .x <- genotypes_of(.x)
  if (is.factor(chromosomes) || is.character(chromosomes)) {
    return(as.character(chromosomes))
  }
  if (!is.numeric(chromosomes)) {
    stop("'chromosomes' should be integer codes of chromosomes")
  }
  chromosome_sidecar_names(.x, chromosomes)
}

#' @rdname gt_lookup_loci
#' @export
gt_restore_loci <- function(.x, names = TRUE, chromosomes = TRUE,
                            check_positions = TRUE) {
  if (!inherits(.x, "gen_tbl")) {
    stop(".x should be a gen_tibble")
  }
  loci <- show_loci(.x)
  if (names && is_int_column(loci$name)) {
    positions <- if (check_positions) loci$position else NULL
    loci$name <- loci_sidecar_names(.x, loci$name, positions)
  }
  if (chromosomes && is_int_column(loci$chromosome)) {
    loci$chromosome <- chromosome_sidecar_factor(.x, loci$chromosome)
  }
  # (a factor of chromosomes is kept as it is by show_loci<-)
  show_loci(.x) <- loci
  .x
}

# Sidecar files of the integer names of loci
#
# When a gen_tibble stores the names of its loci as integers (names_as_int =
# TRUE), the readers write their original names, chromosomes and positions to
# a locus sidecar next to the backing file (foo.bk -> foo.loci.tsv): plain
# tab separated text, with a header, then one line per locus whose first
# field is its integer name (the key, strictly increasing), no quoting, and
# "." for missing values (see src/loci_sidecar.h, and
# data-raw/dev_notes/loci_sidecar_plan.md). The sidecar is found from the
# path of the backing file, so the gen_tibble needs no attribute for it (and
# gt_load(reattach_to = ) finds it too). Sidecars are never edited: functions
# that write a new backing file write a new sidecar for it, with the lines of
# the loci that they keep (copy_loci_sidecar()); integer names keep their
# values, so they stay valid keys.

# the header of a locus sidecar
loci_sidecar_header <- "int_name\tname\tchromosome\tposition"

# the path of a sidecar ("loci" or "chromosomes") of a backing file
sidecar_path_from_bk <- function(bk_file, type = c("loci", "chromosomes")) {
  type <- match.arg(type)
  bigstatsr::sub_bk(bk_file, paste0(".", type, ".tsv"))
}

# the path of the locus sidecar of a gen_tibble (or of its genotypes column)
loci_sidecar_path <- function(x) {
  sidecar_path_from_bk(gt_get_file_names(x)[2], "loci")
}

# the genotypes column of a gen_tibble (or the column itself)
genotypes_of <- function(x) {
  if (inherits(x, "gen_tbl")) {
    x <- x$genotypes
  }
  if (!inherits(x, "vctrs_bigSNP")) {
    stop(".x should be a gen_tibble")
  }
  x
}

# remove the sidecars of a new backing file that are left from older files
# with the same name (which can not belong to the new data), but those that
# are listed in keep
remove_stale_sidecars <- function(bk_file, keep = character(0)) {
  types <- setdiff(c("loci", "chromosomes"), keep)
  for (type in types) {
    path <- sidecar_path_from_bk(bk_file, type)
    if (file.exists(path)) {
      unlink(path)
    }
  }
  invisible(NULL)
}

# the sparse indices of the sidecars read in this session, with the size and
# modification time of each file when it was indexed (a file that changed is
# indexed again)
sidecar_cache <- new.env(parent = emptyenv())

# the sparse index of a sidecar (the key and byte offset of every 1024th
# line), from the cache if the file has not changed
sidecar_get_index <- function(path, header) {
  # (a smaller block is only used by the tests, to cross block edges)
  block <- as.integer(getOption("tidypopgen.sidecar_block", 1024L))
  info <- file.info(path, extra_cols = FALSE)
  cache_key <- normalizePath(path)
  cached <- sidecar_cache[[cache_key]]
  if (!is.null(cached) && identical(cached$size, info$size) &&
        identical(cached$mtime, info$mtime) &&
        identical(cached$block, block)) {
    return(cached$index)
  }
  index <- sidecar_index(path, header, block)
  assign(cache_key,
    list(size = info$size, mtime = info$mtime, block = block, index = index),
    envir = sidecar_cache
  )
  index
}

# integer keys (e.g. the integer names of loci), checking that they are
# positive whole numbers (only checking for whole numbers if they are not
# integers, to avoid creating vectors as long as the loci table)
as_sidecar_keys <- function(keys, what) {
  if (!is.integer(keys)) {
    if (!is.numeric(keys) || any(keys != round(keys), na.rm = TRUE) ||
          any(abs(keys) > .Machine$integer.max, na.rm = TRUE)) {
      stop(what, " should be integers")
    }
    keys <- as.integer(keys)
  }
  if (anyNA(keys) || ((length(keys) > 0) && (min(keys) < 1))) {
    stop(what, " should be positive integers, with no NAs")
  }
  keys
}

# the original names of the loci of x with integer names keys, from the
# locus sidecar of x; with positions, check that the positions in the
# sidecar are the same
loci_sidecar_names <- function(x, keys, positions = NULL) {
  path <- loci_sidecar_path(x)
  if (!file.exists(path)) {
    stop(
      "there is no locus sidecar (", path, ") with the original names of ",
      "the loci of this gen_tibble. It is written next to the backing file ",
      "when a VCF is read with names_as_int = TRUE (but not by older ",
      "versions of tidypopgen, nor yet for other formats); if the files of ",
      "the gen_tibble were moved, move the sidecar with them"
    )
  }
  keys <- as_sidecar_keys(keys, "the names of the loci")
  index <- sidecar_get_index(path, loci_sidecar_header)
  check <- !is.null(positions)
  res <- sidecar_fetch(path, keys,
    index_keys = index$keys, index_offsets = index$offsets,
    fields = 2L, types = "character",
    check_field = if (check) 4L else 0L,
    check_values = if (check) as.integer(positions) else integer(0)
  )
  if (!all(res$found)) {
    not_found <- unique(keys[!res$found])
    stop(
      "the locus sidecar ", path, " has no line for ", length(not_found),
      " of the loci (e.g. the locus with integer name ", not_found[1],
      "): they are not loci of this gen_tibble"
    )
  }
  if (res$n_mismatch > 0) {
    i <- res$first_mismatch
    stop(
      "the locus sidecar ", path, " does not match the loci: ",
      res$n_mismatch, " of ", length(keys), " loci have a different ",
      "position (e.g. the locus with integer name ", keys[i], " is at ",
      "position ", positions[i], ", but at ", res$first_mismatch_value,
      " in the sidecar), so the sidecar probably belongs to other data. If ",
      "you changed the positions of the loci yourself (e.g. after a ",
      "liftover), use check_positions = FALSE"
    )
  }
  res[[1]]
}

# write the locus sidecar of a new backing file (new_bk_file), with the lines
# of the loci of x whose integer names are keys (e.g. the loci kept by
# gt_update_backingfile()). Nothing is written if the names are not integers,
# or if x has no sidecar (e.g. data read with older versions of tidypopgen,
# whose lookups give a clear error); if the sidecar of x lacks some of the
# loci, the new sidecar is removed, with a warning. It returns whether a
# sidecar was written
copy_loci_sidecar <- function(x, new_bk_file, keys) {
  new_path <- sidecar_path_from_bk(new_bk_file, "loci")
  # (a sidecar left by older files with the same name would not match)
  if (file.exists(new_path)) {
    unlink(new_path)
  }
  path <- loci_sidecar_path(x)
  if (!is.numeric(keys) || !file.exists(path)) {
    return(invisible(FALSE))
  }
  keys <- as_sidecar_keys(keys, "the names of the loci")
  res <- sidecar_subset(path, new_path, loci_sidecar_header, keys)
  if (res$n_written < res$n_keys) {
    unlink(new_path)
    warning(
      "the locus sidecar ", path, " does not have all the loci, so the ",
      "original names of the loci of the new gen_tibble are not kept"
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# the chromosome sidecar of x has the names of its integer chromosomes
# (chromosomes_as_int): foo.bk -> foo.chromosomes.tsv, with a line for each
# chromosome (its code and name), in increasing order of code

# the header of a chromosome sidecar
chromosomes_sidecar_header <- "int_chromosome\tchromosome"

# the path of the chromosome sidecar of a gen_tibble (or of its genotypes
# column)
chromosomes_sidecar_path <- function(x) {
  sidecar_path_from_bk(gt_get_file_names(x)[2], "chromosomes")
}

# whether a column of the loci table holds integers (integer names or
# chromosome codes), rather than names (a character vector or a factor)
is_int_column <- function(x) {
  is.numeric(x) && !is.factor(x)
}

# stop if the chromosome sidecar at path does not exist
check_chromosome_sidecar <- function(path) {
  if (!file.exists(path)) {
    stop(
      "there is no chromosome sidecar (", path, ") with the names of the ",
      "chromosomes of this gen_tibble. It is written next to the backing ",
      "file when a VCF is read with chromosomes_as_int = TRUE; if the files ",
      "of the gen_tibble were moved, move the sidecar with them"
    )
  }
  invisible(NULL)
}

# the names of the chromosome codes of x (NA for NA), from the chromosome
# sidecar of x; each distinct code is looked up once
chromosome_sidecar_names <- function(x, codes) {
  path <- chromosomes_sidecar_path(x)
  check_chromosome_sidecar(path)
  values <- unique(codes)
  values <- values[!is.na(values)]
  keys <- as_sidecar_keys(values, "the chromosome codes")
  index <- sidecar_get_index(path, chromosomes_sidecar_header)
  res <- sidecar_fetch(path, keys,
    index_keys = index$keys, index_offsets = index$offsets,
    fields = 2L, types = "character",
    check_field = 0L, check_values = integer(0)
  )
  if (!all(res$found)) {
    stop(
      "the chromosome sidecar ", path, " has no line for ",
      sum(!res$found), " of the chromosome codes (e.g. ",
      keys[!res$found][1], "): they are not chromosomes of this gen_tibble"
    )
  }
  res[[1]][match(codes, values)]
}

# the chromosome codes of x as a factor of their names (as gen_tibble()
# gives without chromosomes_as_int: the levels are the names in the
# chromosome sidecar, sorted as as.factor() sorts them), reading the whole
# sidecar (one line per chromosome) and creating no string per locus
chromosome_sidecar_factor <- function(x, codes) {
  path <- chromosomes_sidecar_path(x)
  check_chromosome_sidecar(path)
  index <- sidecar_get_index(path, chromosomes_sidecar_header)
  sidecar <- sidecar_read(path, chromosomes_sidecar_header,
    fields = c(1L, 2L), types = c("integer", "character"), n = index$n
  )
  chromosome_index <- match(codes, sidecar[[1]])
  not_found <- is.na(chromosome_index) & !is.na(codes)
  if (any(not_found)) {
    stop(
      "the chromosome sidecar ", path, " has no line for some of the ",
      "chromosome codes (e.g. ", codes[not_found][1], "): they are not ",
      "chromosomes of this gen_tibble"
    )
  }
  chromosome_index_to_factor(chromosome_index, sidecar[[2]])
}

# write the chromosome sidecar of a new backing file (new_bk_file), with the
# lines of the chromosome codes of x in codes (see copy_loci_sidecar())
copy_chromosomes_sidecar <- function(x, new_bk_file, codes) {
  new_path <- sidecar_path_from_bk(new_bk_file, "chromosomes")
  if (file.exists(new_path)) {
    unlink(new_path)
  }
  path <- chromosomes_sidecar_path(x)
  if (!is_int_column(codes) || !file.exists(path)) {
    return(invisible(FALSE))
  }
  keys <- as_sidecar_keys(unique(codes[!is.na(codes)]), "the chromosome codes")
  res <- sidecar_subset(path, new_path, chromosomes_sidecar_header, keys)
  if (res$n_written < res$n_keys) {
    unlink(new_path)
    warning(
      "the chromosome sidecar ", path, " does not have all the chromosomes, ",
      "so the names of the chromosomes of the new gen_tibble are not kept"
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# write the sidecars of a new backing file (new_bk_file) whose loci table is
# loci, with the loci (and chromosomes) of x that it keeps
copy_sidecars <- function(x, new_bk_file, loci) {
  copy_loci_sidecar(x, new_bk_file, loci$name)
  copy_chromosomes_sidecar(x, new_bk_file, loci$chromosome)
  invisible(NULL)
}

# the sidecars of a gen_tibble that exist, for the messages of gt_save()
sidecar_files <- function(x) {
  loci <- show_loci(x)
  paths <- character(0)
  if (is_int_column(loci$name)) {
    paths[["loci"]] <- loci_sidecar_path(x)
  }
  if (is_int_column(loci$chromosome)) {
    paths[["chromosomes"]] <- chromosomes_sidecar_path(x)
  }
  paths[file.exists(paths)]
}

# message the sidecars of a gen_tibble (after the files of its FBM)
message_sidecar_files <- function(x) {
  paths <- sidecar_files(x)
  if ("loci" %in% names(paths)) {
    message("with original names of loci in: ", paths[["loci"]])
  }
  if ("chromosomes" %in% names(paths)) {
    message("with names of chromosomes in: ", paths[["chromosomes"]])
  }
  invisible(NULL)
}
