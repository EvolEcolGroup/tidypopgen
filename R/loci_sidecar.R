#' Look up the original names of loci stored as integers
#'
#' With `names_as_int = TRUE`, [gen_tibble()] stores the names of the loci as
#' integers, to save memory with large datasets, and writes their original
#' names, chromosomes and positions to a sidecar file next to the backing file
#' (see Details). `gt_lookup_loci()` gives back the original names of some
#' loci (e.g. the hits of a selection scan), without changing the
#' `gen_tibble`; `gt_restore_loci()` puts the original names back into the
#' loci table of a `gen_tibble` (e.g. once it has been filtered down to fewer
#' loci).
#'
#' @details The sidecar of the backing file `foo.bk` is `foo.loci.tsv`: a tab
#'   separated text file with a header (`int_name`, `name`, `chromosome`,
#'   `position`) and one line per locus, in increasing order of integer name.
#'   It has no quotes, and "." is a missing value, so it can be read with any
#'   tool (e.g. `data.table::fread(path, quote = "")` or
#'   `readr::read_tsv(path, quote = "")`), and joined to [show_loci()] by
#'   `name = int_name`. The names are those that `names_as_int = FALSE` would
#'   give (with a missing ID replaced by the chromosome and position).
#'
#'   The sidecar is written by [gen_tibble()] when it reads a VCF file with
#'   `parser = "cpp"`. Functions that write a new backing file
#'   ([gt_update_backingfile()], [gt_order_loci()], [rbind()][rbind.gen_tbl])
#'   write a new sidecar for it, with the loci that they keep (integer names
#'   keep their values, so they still find their line). Keep the sidecar with
#'   the other files of the `gen_tibble` (see [gt_save()]): without it, the
#'   original names are lost.
#'
#'   Only the file name ties a sidecar to its data, so the position of each
#'   locus in the sidecar is compared with its position in the loci table: a
#'   different position means that the sidecar belongs to other data (e.g.
#'   files copied by hand from different runs), and is an error. If you
#'   changed the positions of the loci yourself (e.g. after a liftover to a
#'   new assembly), the names in the sidecar are still right: use
#'   `check_positions = FALSE`.
#'
#'   Lookups are fast: the first lookup in a session reads through the
#'   sidecar once to index it (about a second per GB), and later lookups of a
#'   few loci only read the parts of the file that hold them.
#'
#' @param .x a [`gen_tibble`] (or, for `gt_lookup_loci()`, its `genotypes`
#'   column)
#' @param loci a subset of the loci table of `.x`: rows of [show_loci()], in
#'   any order, possibly with extra columns (e.g. the statistics of a scan);
#'   or row numbers of `show_loci(.x)`, as given by [which()].
#' @param check_positions boolean: whether to check that the positions of the
#'   loci in the sidecar are the same as in `loci` (or in the loci table of
#'   `.x`, for `gt_restore_loci()`). For `gt_lookup_loci()`, they are only
#'   checked if `loci` has a `position` column. Default is TRUE.
#' @returns `gt_lookup_loci()` returns `loci` (as a table), with the integer
#'   names in its `name` column replaced by the original names; the other
#'   columns, and the rows and their order, are unchanged. A `name` column
#'   that already holds names is left as it is, so the same code works on data
#'   read with or without `names_as_int`.
#'
#'   `gt_restore_loci()` returns `.x`, with the original names in its loci
#'   table (as if it had been read with `names_as_int = FALSE`). It does not
#'   save it: use [gt_save()].
#' @rdname gt_lookup_loci
#' @export
#' @examples
#' vcf_path <- system.file(
#'   "extdata/anolis/punctatus_t70_s10_n46_filtered.recode.vcf.gz",
#'   package = "tidypopgen"
#' )
#' anole_gt <- gen_tibble(vcf_path,
#'   backingfile = tempfile("anolis_"),
#'   names_as_int = TRUE, quiet = TRUE
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
  if (!"name" %in% names(loci)) {
    stop("'loci' should have a 'name' column, as the loci table of .x")
  }
  if (is.numeric(loci$name)) {
    positions <- NULL
    if (check_positions && ("position" %in% names(loci))) {
      positions <- loci$position
    }
    loci$name <- loci_sidecar_names(.x, loci$name, positions)
  }
  loci
}

#' @rdname gt_lookup_loci
#' @export
gt_restore_loci <- function(.x, check_positions = TRUE) {
  if (!inherits(.x, "gen_tbl")) {
    stop(".x should be a gen_tibble")
  }
  loci <- show_loci(.x)
  if (is.numeric(loci$name)) {
    positions <- if (check_positions) loci$position else NULL
    show_loci(.x)$name <- loci_sidecar_names(.x, loci$name, positions)
  }
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

# the sidecars of a gen_tibble that exist, for the messages of gt_save()
sidecar_files <- function(x) {
  path <- loci_sidecar_path(x)
  if (is.numeric(show_loci(x)$name) && file.exists(path)) {
    return(path)
  }
  character(0)
}

# message the sidecars of a gen_tibble (after the files of its FBM)
message_sidecar_files <- function(x) {
  for (path in sidecar_files(x)) {
    message("with original names of loci in: ", path)
  }
  invisible(NULL)
}
