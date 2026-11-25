# time the steps of gen_tibble_packedancestry for the new version
args <- commandArgs(TRUE); lib <- args[1]; file <- args[2]; n_cores <- as.integer(args[3]); bk <- args[4]
suppressMessages(library(tidypopgen, lib.loc = lib))
ns <- asNamespace("tidypopgen")
t_all <- system.time(gt <- gen_tibble(file, backingfile = bk, n_cores = n_cores, quiet = TRUE))
map_file <- sub("\\.t?geno$", ".snp", file)
t_snp <- system.time(lt <- ns$read_loci_file(map_file, character_cols = c(1, 5, 6)))
names(lt) <- c("name", "chromosome", "genetic_dist", "position", "allele_ref", "allele_alt")
t_val <- system.time(ns$validate_loci(lt, check_alphabet = TRUE, harmonise_loci = TRUE, check_duplicates = TRUE, allow_duplicates = FALSE, valid_alleles = c("A","T","C","G"), missing_alleles = c("0", ".")))
cat(sprintf("STEPS file=%s cores=%d total=%.2f (user %.2f sys %.2f) | read .snp %.2f | validate_loci %.2f\n",
  basename(file), n_cores, t_all[["elapsed"]], t_all[["user.self"]], t_all[["sys.self"]], t_snp[["elapsed"]], t_val[["elapsed"]]))
