# Time the phases of gt_merge_vcfs(): a copy of its body (and of
# write_vcfs_gen_tibble()) with timers around each step, run with the
# internals of an installed tidypopgen.
#   Rscript merge_phases.R <lib> <scenario> <n_cores> <backingfile> <out_tsv>
#     [real]
# scenario: pops (ref_mismatch = "swap", intersection), sgdp / many (union,
# missing_to_ref = TRUE). With "real" as the last argument, it calls
# gt_merge_vcfs() itself (no phases), to check that the copy is faithful.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
scenario <- args[2]
n_cores <- as.integer(args[3])
bk <- args[4]
out_tsv <- args[5]
real <- length(args) >= 6 && args[6] == "real"
# lib is a library with tidypopgen installed, or "loadall:<path>" to load a
# checkout compiled at -O2 (pkgbuild::compile_dll(debug = FALSE))
if (startsWith(lib, "loadall:")) {
  pkgload::load_all(sub("^loadall:", "", lib), compile = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(tidypopgen, lib.loc = lib))
}
data_dir <- Sys.getenv(
  "MERGE_DATA", "/media/andrea/annex/tmp/claude/tidypopgen/merge_vcf/data"
)
files <- sort(list.files(file.path(data_dir, scenario),
  pattern = "\\.vcf\\.gz$", full.names = TRUE
))
opts <- switch(scenario,
  pops = list(variants = "intersection", missing_to_ref = FALSE,
              ref_mismatch = "swap"),
  list(variants = "union", missing_to_ref = TRUE, ref_mismatch = "error")
)
names_as_int <- Sys.getenv("MERGE_NAMES_AS_INT", "TRUE") == "TRUE"

ns <- asNamespace("tidypopgen")
for (fn in c(
  "check_vcf_files_args", "default_vcfs_backingfile", "vcfs_sidecar_paths",
  "merge_sites_new", "scan_vcf_file", "merge_vcfs_sample_names",
  "merge_sites_add_file", "sidecar_index", "loci_sidecar_header",
  "sidecar_block", "sidecar_gather", "merge_sites_result",
  "stop_ref_mismatch", "report_vcf_merge", "chromosome_index_to_factor",
  "validate_loci", "merge_sites_file_cols", "validate_indiv_meta",
  "vcf_merge_genotypes", "new_vctrs_bigsnp", "has_genotypes_without_alleles",
  "gt_save_light", "remove_fbm_files"
)) {
  assign(fn, get(fn, envir = ns))
}
# features of the prototypes (mirrored below)
has_cache <- exists("bk_merge_genotypes", envir = ns) &&
  getOption("tidypopgen.merge_cache", TRUE)
if (has_cache) bk_merge_genotypes <- get("bk_merge_genotypes", envir = ns)
# (alleles as raw vectors, and the one-pass gather of new loci)
has_raw <- "alleles_as_raw" %in% names(formals(get("vcf_to_bk", envir = ns)))
cat("has_cache:", has_cache, " has_raw:", has_raw, "\n")

# accumulated elapsed time of each phase
phase_times <- list()
tick <- function() proc.time()[["elapsed"]]
add_time <- function(phase, t0) {
  phase_times[[phase]] <<- (if (is.null(phase_times[[phase]])) 0 else
    phase_times[[phase]]) + tick() - t0
}

t_start <- tick()
if (real) {
  x_gt <- gt_merge_vcfs(files,
    variants = opts$variants, missing_to_ref = opts$missing_to_ref,
    ref_mismatch = opts$ref_mismatch, backingfile = bk, n_cores = n_cores,
    names_as_int = names_as_int, quiet = TRUE
  )
  add_time("total_real", t_start)
} else {
  quiet <- TRUE
  variants <- opts$variants
  missing_to_ref <- opts$missing_to_ref
  alt_mismatch <- "drop"
  ref_mismatch <- opts$ref_mismatch
  force_samples <- FALSE
  valid_alleles <- c("A", "T", "C", "G")
  missing_alleles <- c("0", ".")
  allow_duplicates <- FALSE
  chromosomes_as_int <- FALSE
  phased <- FALSE
  files <- check_vcf_files_args(files, valid_alleles, n_cores)
  backingfile <- default_vcfs_backingfile(bk, files, "_merged")
  n_files <- length(files)
  sidecars <- vcfs_sidecar_paths(backingfile)
  n_sites <- 0
  merged_sites <- merge_sites_new(
    split_alt = alt_mismatch == "split",
    ref_rule = match(ref_mismatch, c("error", "drop", "swap")) - 1L,
    missing_alleles = union(c(missing_alleles, "0"), ".")
  )
  sample_names <- vector("list", n_files)
  file_ploidy <- vector("list", n_files)
  site_names <- vector("list", n_files)
  chromosome_names <- character(0)
  cache_files <- if (has_cache) {
    paste0(backingfile, "_cache", seq_len(n_files), ".bk")
  } else {
    character(0)
  }
  cache_missing <- integer(n_files)
  cache_swap_ok <- logical(n_files)
  for (f in seq_len(n_files)) {
    t0 <- tick()
    scan <- if (has_cache) {
      if (has_raw) {
        scan_vcf_file(files[f], names_as_int, n_cores, quiet,
          loci_file = if (names_as_int) sidecars$scan else "",
          bk_file = cache_files[f], phased = phased, alleles_as_raw = TRUE
        )
      } else {
        scan_vcf_file(files[f], names_as_int, n_cores, quiet,
          loci_file = if (names_as_int) sidecars$scan else "",
          bk_file = cache_files[f], phased = phased
        )
      }
    } else {
      scan_vcf_file(files[f], names_as_int, n_cores, quiet,
        loci_file = if (names_as_int) sidecars$scan else ""
      )
    }
    add_time("p1_scan_cpp", t0)
    t0 <- tick()
    sample_names[[f]] <- scan$sample_names
    file_ploidy[[f]] <- scan$ploidy
    if (has_cache) {
      cache_missing[f] <- max(scan$ploidy) + 1L
      cache_swap_ok[f] <- is.na(scan$min_alleles) ||
        (scan$min_alleles == scan$max_alleles &&
          all(scan$ploidy == scan$max_alleles))
    }
    merge_vcfs_sample_names(sample_names[seq_len(f)], force_samples)
    chromosome_names <- vctrs::vec_unique(
      c(chromosome_names, scan$chromosome_names)
    )
    chr_codes <- match(scan$chromosome_names, chromosome_names)[
      scan$loci_tbl$chromosome
    ]
    add_time("p1_r_glue", t0)
    t0 <- tick()
    new_loci <- merge_sites_add_file(merged_sites,
      chromosome = chr_codes,
      position = scan$loci_tbl$physical.pos,
      allele_ref = scan$loci_tbl$allele2,
      allele_alt = scan$loci_tbl$allele1
    )
    add_time("p1_match_cpp", t0)
    t0 <- tick()
    if (!names_as_int) {
      site_names[[f]] <- scan$loci_tbl$marker.ID[new_loci]
    } else if (length(new_loci) > 0) {
      scan_index <- if (has_raw) {
        list(keys = integer(0), offsets = numeric(0))
      } else {
        sidecar_index(sidecars$scan, loci_sidecar_header, sidecar_block())
      }
      sidecar_gather(sidecars$scan, loci_sidecar_header,
        index_keys = scan_index$keys, index_offsets = scan_index$offsets,
        keys = new_loci, path_out = sidecars$sites,
        first_new_key = n_sites + 1, append = file.exists(sidecars$sites)
      )
    }
    unlink(sidecars$scan)
    n_sites <- n_sites + length(new_loci)
    rm(scan, new_loci, chr_codes)
    add_time("p1_sidecar", t0)
  }
  t0 <- tick()
  sites <- merge_sites_result(merged_sites)
  add_time("mid_result", t0)
  t0 <- tick()
  if (sites$ref_mismatch$n > 0) {
    stop_ref_mismatch(sites, chromosome_names, files)
  }
  sample_id <- merge_vcfs_sample_names(sample_names, force_samples)
  n_indiv <- lengths(sample_names)
  keep <- sites$status == 0L
  if (variants == "intersection") {
    keep <- keep & sites$n_files == n_files
  }
  keep <- which(keep)
  keep <- keep[order(
    sites$chromosome[keep], sites$position[keep],
    method = "radix"
  )]
  chromosome_used <- rep(FALSE, length(chromosome_names))
  chromosome_used[sites$chromosome[keep]] <- TRUE
  chromosome_new_index <- cumsum(chromosome_used)
  used_names <- chromosome_names[chromosome_used]
  loci_chromosome <- chromosome_index_to_factor(
    chromosome_new_index[sites$chromosome[keep]],
    used_names
  )
  loci_position <- sites$position[keep]
  add_time("mid_keep_order", t0)
  t0 <- tick()
  if (names_as_int) {
    loci_name <- seq_along(keep)
    sites_index <- sidecar_index(
      sidecars$sites, loci_sidecar_header, sidecar_block()
    )
    sidecar_gather(sidecars$sites, loci_sidecar_header,
      index_keys = sites_index$keys, index_offsets = sites_index$offsets,
      keys = keep, path_out = sidecars$loci,
      first_new_key = 1, append = FALSE
    )
    unlink(sidecars$sites)
  } else {
    loci_name <- unlist(site_names)[keep]
    rm(site_names)
  }
  add_time("mid_sidecar", t0)
  t0 <- tick()
  loci <- tibble::tibble(
    name = loci_name,
    chromosome = loci_chromosome,
    position = loci_position,
    genetic_dist = 0,
    allele_ref = sites$allele_ref[keep],
    allele_alt = sites$allele_alt[keep]
  )
  loci <- validate_loci(loci,
    check_alphabet = TRUE,
    harmonise_loci = TRUE,
    check_duplicates = TRUE,
    allow_duplicates = allow_duplicates,
    valid_alleles = valid_alleles,
    missing_alleles = missing_alleles,
    chromosomes_as_int = chromosomes_as_int
  )
  site_col <- rep(-1L, length(sites$status))
  site_col[keep] <- seq_along(keep) - 1L
  rm(sites)
  add_time("mid_validate_loci", t0)

  # write_vcfs_gen_tibble()
  t0 <- tick()
  fill_missing <- variants == "union" && !missing_to_ref
  ploidy <- unlist(file_ploidy)
  row_offset <- c(0L, cumsum(n_indiv))[seq_len(n_files)]
  indiv_meta <- validate_indiv_meta(data.frame(id = sample_id))
  max_ploidy <- max(ploidy)
  missing_value <- max_ploidy + 1
  code256 <- rep(NA_real_, 256)
  code256[1:(max_ploidy + 1)] <- seq(0, max_ploidy)
  file_backed_matrix <- bigstatsr::FBM.code256(
    nrow = length(sample_id),
    ncol = nrow(loci),
    code = code256,
    backingfile = backingfile,
    init = if (fill_missing) missing_value else NULL
  )
  bk_file <- file_backed_matrix$backingfile
  add_time("p2_fbm_create", t0)
  for (f in seq_along(files)) {
    t0 <- tick()
    this_file <- merge_sites_file_cols(merged_sites, f, site_col)
    add_time("p2_file_cols", t0)
    t0 <- tick()
    if (has_cache && (!any(this_file$swap) || cache_swap_ok[f])) {
      bk_merge_genotypes(cache_files[f],
        BM = file_backed_matrix,
        row_offset = row_offset[f],
        row = integer(0),
        col = this_file$col,
        swap = this_file$swap,
        ploidy = file_ploidy[[f]],
        cache_missing = cache_missing[f],
        missing_value = missing_value,
        n_threads = n_cores
      )
    } else {
      vcf_merge_genotypes(files[f],
        BM = file_backed_matrix,
        row_offset = row_offset[f],
        row = integer(0),
        col = this_file$col,
        swap = this_file$swap,
        ploidy = file_ploidy[[f]],
        missing_value = missing_value,
        n_threads = n_cores,
        phased_lut = integer(0),
        sample_names = character(0)
      )
    }
    if (has_cache) unlink(cache_files[f])
    add_time("p2_genotypes_cpp", t0)
  }
  t0 <- tick()
  fbm_path <- bigstatsr::sub_bk(bk_file, ".rds")
  saveRDS(file_backed_matrix, fbm_path)
  indiv_meta$genotypes <- new_vctrs_bigsnp(
    fbm_obj = file_backed_matrix,
    fbm_file = fbm_path,
    loci = loci,
    indiv_id = indiv_meta$id,
    ploidy = max_ploidy,
    fbm_ploidy = ploidy
  )
  x_gt <- tibble::new_tibble(indiv_meta, class = "gen_tbl")
  rm(file_backed_matrix)
  add_time("end_new_gt", t0)
  t0 <- tick()
  stopifnot(!has_genotypes_without_alleles(x_gt))
  add_time("end_check_alleles", t0)
  t0 <- tick()
  gt_save_light(x_gt, quiet = quiet)
  add_time("end_save_gt", t0)
}
add_time("total", t_start)
# the samples and loci, to compare the outputs of versions
saveRDS(
  list(id = x_gt$id, loci = show_loci(x_gt)[, -1], ploidy = indiv_ploidy(x_gt)),
  paste0(bk, "_meta.rds")
)
res <- data.frame(
  scenario = scenario, n_cores = n_cores, rep = Sys.getenv("MERGE_REP", "0"),
  phase = names(phase_times),
  seconds = round(unlist(phase_times), 3), n_indiv = nrow(x_gt),
  n_loci = count_loci(x_gt), row.names = NULL
)
write.table(res, out_tsv,
  sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = !file.exists(out_tsv), append = TRUE
)
print(res)
