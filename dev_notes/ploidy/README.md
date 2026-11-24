# Arbitrary ploidy: roadmap and checklist

The aim is that every analysis in tidypopgen either works on haploid,
polyploid and mixed-ploidy data, or says clearly that it does not. This is
the checklist; the reasoning (estimators, validation, design choices) is in
`ploidy_roadmap.typ` and its PDF.

A box is ticked when the work is on `dev`. Items that are written but not
merged name their branch. Each branch is one small step, off `dev`;
`ploidy_guards` is the base of the other code branches and merges first.

`audit.R` calls every analysis on diploid, haploid, tetraploid, hexaploid and
mixed (2x + 4x) data and writes `audit.tsv` (`audit_v1.0.0.tsv` is the
baseline). It is the measure of progress: rerun it as branches merge.

## Already in v1.0.0

- [x] Dosage storage with a ploidy for each individual; `show_ploidy()`,
      `indiv_ploidy()`; the label is refreshed by `filter()`
- [x] `gen_tibble()` reads polyploid and mixed-ploidy VCF files and matrices
- [x] `loci_missingness()`, `indiv_missingness()`
- [x] `indiv_het_obs()`, `pop_het_obs()`
- [x] `loci_pi()`
- [x] `pop_tajimas_d()`, with windows and jackknife
- [x] `pop_sfs()`, `pairwise_pop_sfs()`, `nwise_pop_sfs()`
- [x] `gt_as_genind()`
- [x] "diploid only" errors in `gt_as_hierfstat()`, `gt_as_geno_lea()`,
      `gt_as_plink()` and most diploid statistics

## Phase 0: foundations

- [ ] One classifier and one guard, `.ploidy_class()` and `.check_ploidy()`;
      errors name the function and the ploidy (branch `ploidy_guards`)
- [ ] Test helpers: simulation for any ploidy, reference counts, the merging
      invariant (branch `ploidy_guards`)
- [ ] PCA wrappers, imputation, pcadapt and ROH stop on unsupported ploidy
      instead of treating dosages as diploid (branch
      `ploidy_guard_unguarded`)
- [ ] The audit as a snapshot test
- [ ] Rerun the audit on `dev` (the tables here predate the rebase onto it)

## Phase 1: allele frequencies and gene diversity

- [ ] Grouped `loci_alt_freq()` and `loci_maf()` (branch `ploidy_loci_freq`)
- [ ] `qc_report_loci()`: MAF and missingness for any ploidy; fixes halved
      haploid frequencies (branch `ploidy_qc_report_loci`)
- [ ] `pop_het_exp()`, `pop_gene_div()` (branch `ploidy_pop_het_exp`)
- [ ] `pop_het_exp(include_global = TRUE)` (with `pop_global_stats()`)

## Phase 2: differentiation between populations

- [ ] Hudson's Fst: `pairwise_pop_fst()`, windows, jackknife, PBS (branch
      `ploidy_hudson_fst`)
- [ ] `pairwise_allele_sharing()`
- [ ] `pop_fst()` (WG17)
- [ ] `pop_fis(method = "WG17")`, `indiv_inbreeding()`
- [ ] `pairwise_grm()` (WG17)
- [ ] Nei87: `pairwise_pop_fst()`, `pop_fis()`
- [ ] `pop_global_stats()`
- [ ] `method = "rho"` for polyploids (WC84 stays diploid-only)
- [ ] `gt_extract_f2()` checked against ADMIXTOOLS 2

## Phase 3: relationships between individuals

- [ ] `pairwise_ibs()` on the frequency scale
- [ ] VanRaden GRM for autopolyploids
- [ ] `qc_report_indiv()` and `filter_high_relatedness()` with allele sharing
      kinship (KING stays diploid-only)

## Phase 4: multivariate analyses and imputation

- [ ] Re-read this phase against `gt_pca()` (new in 1.1.0)
- [ ] PCA for a single ploidy
- [ ] PCA for mixed ploidy (design and benchmark first)
- [ ] Imputation codes for any ploidy; `gt_impute_simple()`
- [ ] `gt_dapc()`, `gt_pcadapt()`, `loci_ld_clump()`, `loci_ld_prune()`
- [ ] `gt_snmf()` takes the ploidy from the data

## Phase 5: input, output and combining data

- [ ] `gt_as_vcf()`
- [ ] `rbind()`
- [ ] `gt_as_genlight()`, `gt_from_genlight()`
- [ ] STRUCTURE export (to be decided)

## Real data

Everything so far is tested on simulated dosages. A real dataset is needed to
find what simulations miss (allele dosage errors, uneven missingness between
ploidies, VCF conventions of polyploid callers), for the vignette, and for
benchmarks. It can be chosen at any time, and the sooner the better.

- [ ] Choose a real polyploid dataset. What to look for: a public VCF (or
      reads and a reference) with a licence that allows a subset in the
      package or its articles; called with a polyploid-aware caller, with
      dosages rather than diploidised genotypes; several populations with
      enough individuals each; ideally diploids and polyploids of the same
      species, to test mixed ploidy; and published statistics (Fst, diversity,
      PCA) to compare against
- [ ] First candidate, for tests, examples and the vignette: the *Arabidopsis
      arenosa* tutorial VCF of Bohutinska et al. 2023 (Methods Mol Biol
      2545:297-324), in `github.com/mbohutinska/PolyplChapter` (`PartB`). 40
      individuals in 5 populations (3 diploid, 2 autotetraploid), 69,245
      biallelic SNPs on one scaffold, GATK calls with real dosages (file
      inspected, 25 Nov 2026); the repository has per-population pi, Tajima's
      D, Fst, rho and dxy to compare with. The repository has no licence:
      ask the authors before putting a subset in the package
- [ ] Second candidate, for an article and benchmarks: the *Cardamine* 2x/4x
      (and 3x hybrids?) ddRAD VCF on Dryad, doi:10.5061/dryad.8931zcs69 (CC0),
      220 individuals, 456,239 biallelic SNPs, 888 MB, GATK with a ploidy for
      each sample. The file itself has not been opened yet (dosage encoding,
      whether the triploids are in it, ploidy of each sample)
- Looked at and set aside (25 Nov 2026): *A. arenosa* range-wide (Monnahan et
  al. 2019 and later) and switchgrass (reads only, no VCF found); *Cardamine
  amara* (pool-seq); *Cochlearia* (no VCF deposited); potato panels (one
  ploidy, no populations). Fallbacks: *Miscanthus sacchariflorus* 2x/4x (CC
  BY, format not checked) and *Luzula* 2x/4x/6x (1,002 samples, CC0)
- [ ] A small subset in `inst/extdata` (or built in `data-raw`) for examples
      and tests, and the full dataset for an article and for benchmarks
- [ ] Run the audit and the published statistics on it

## Phase 6: validation and documentation

- [ ] Fixtures from other packages in `data-raw/validate`: StAMPP, AGHmatrix,
      GenoDive, adegenet/hierfstat through `gt_as_genind()`
- [ ] Vignette on polyploid and mixed-ploidy data
- [ ] A "Ploidy" section in each help page, and the support table as an
      article

## Later (they change the data model)

- [ ] Ploidy that varies along the genome: sex chromosomes, haplodiploids,
      organelles
- [ ] Uncertain dosages (updog, polyRAD)

## Open decisions

- `loci_hwe()`: stay diploid-only? (proposed: yes, as the inheritance mode is
  not recorded)
- `pairwise_pop_fst()`: default to Hudson for data that are not diploid?
- `rho` or a polyploid WC84?
- NEWS: the branches add a "development version" section above 1.1.0
