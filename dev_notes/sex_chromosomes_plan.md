# Plan: sex chromosomes and mtDNA

Teach tidypopgen which loci are autosomal, on a sex chromosome (X, Y, Z, W
and the pseudo-autosomal regions) or on the mitochondrial genome, and
which individuals are male or female. With that it can:

- store hemizygous genotypes the same way whatever the input file;
- run the analyses that assume diploid autosomes on autosomes only, as
  PLINK does;
- check the recorded sex against the genotypes;
- write sex chromosomes back correctly to PLINK and VCF files.

Status: phase 1 implemented on 26 Nov 2026 (laptopX1), on branch
`sex_chromosomes` (from `dev` at 4374e70d; pushed on 27 Nov, no PR until
Andrea's review), one
commit per step: chromosome types (`b1fb3203`), VCF readers (`9526808e`),
`indiv_sex_check()` (`69bd0209`), autosome-only analyses (`c56ad0bd`),
writers (`815611fa`), article and cheatsheet (`a7a81f4d`). Phase 2 is
outlined only. See "As built" for the differences from the design below.

## Why

tidypopgen knows ploidy per individual only: `attr(genotypes, "ploidy")`
(2; 0 = mixed; -2 = pseudohaploid) and `attr(genotypes, "fbm_ploidy")` (one
value per FBM row). On sex chromosomes and mtDNA, ploidy depends on the
individual *and* the locus. In XY systems males are hemizygous on X and Y
and females have no Y; in ZW systems (birds, Lepidoptera) it is the reverse;
mtDNA is haploid in everyone; the PARs behave like autosomes. On `dev`
(4374e70d):

- **Inconsistent coding.** A male X genotype read from `.bed`, `.ped` or
  EIGENSOFT files is homozygous 0/2, as PLINK stores it. The VCF readers
  store a haploid call as its ALT count, 0/1 (`count_alt_alleles` in
  `src/vcf_parser_fast.cpp`; the vcfR route, `gt_merge_vcfs()` and
  `gt_concat_vcfs()` do the same), so in a diploid individual it looks like
  a heterozygote.
- **Ploidy guessed wrong.** It is guessed from the first VCF line. If that
  line is on X, males silently get ploidy 1 and the tibble becomes mixed
  ploidy; if all samples are haploid there, reading fails. Phased reading
  rejects any haploid call. `gt_pseudohaploid()` tests the first 10,000
  loci, so in an X-first file males would be flagged as pseudohaploid.
- **All loci used everywhere.** No function restricts to autosomes, and
  `gt_extract_f2()` hard-codes `auto_only <- FALSE`. X in males inflates F,
  inbreeding and ROH and fails HWE; it adds sex structure to PCA, KING and
  GRM; it biases frequencies and sample sizes.
- **Chromosome codes clash with PLINK.** `cast_chromosome_to_int()` codes
  non-numeric names after the largest number, in level order, so for a
  factor MT=23 and X=24. `gt_as_plink(chromosomes_as_int = TRUE)`, which
  `gt_admixture()` uses, therefore writes MT as 23, which PLINK reads as X.
  PLINK's own numeric 23-26 codes are taken as autosomes.
- **Sex is only partly read.** `sex` is a factor (female/male) read from
  `.fam`/`.ped`. The packedancestrymap reader reads the `.ind` sex column
  and drops it; VCFs have none.

## What PLINK and EIGENSOFT do (the references)

- **Codes.** With n autosomes, PLINK codes X as n+1, Y as n+2, XY (the PAR)
  as n+3 and MT as n+4: 23-26 for humans (`--chr-set`, `--dog` and so on).
  PLINK 2 prefers PAR1 and PAR2 to XY. EIGENSOFT (CONVERTF README): "X
  chromosome is encoded as 23. Also, Y is encoded as 24, mtDNA is encoded as
  90, and XY is encoded as 91."
- **Haploid calls.** Stored as homozygous (0/2). Heterozygous haploid calls
  and female Y calls are "treated as missing by all analysis commands" in
  PLINK 1.9; `--set-hh-missing` erases them.
- **`--split-x`.** The two numbers are the last position of PAR1 and the
  first of PAR2: b36/hg18 2709521 and 154584237; b37/hg19 2699520 and
  154931044; b38/hg38 2781479 and 155701383.
- **`--check-sex`.** X inbreeding F, "based on the same F coefficient as
  --het/--ibc": F = (O - E) / (N - E), with E summing 1 - 2pq without a
  small-sample correction. Females are called below 0.2 and males above 0.8
  (1.9). `ycount` also uses the number of non-missing Y genotypes. The PARs
  should be split first.
- **Autosome-only commands.**
  - KING (2.0), `--genome` and `--het` (1.9), and Weir-Cockerham `--fst`
    use autosomes only.
  - admixtools' `extract_f2()` does too by default (`auto_only = TRUE`).
  - `--freq` (2.0) counts one allele per male on X.
  - HWE on X is a separate test in PLINK 2 (`.hardy.x`, Graffelman & Weir
    2016).

## Decisions (Andrea, 26 Nov 2026)

0. **Behaviour changes.** Results change for data with loci recognised as
   X, Y or MT (VCF haploid calls stored as 0/2, autosome-only analyses).
   This is accepted without deprecated functions, as the new behaviour is
   the correct one; declaring a chromosome as an autosome brings the old
   behaviour back (Andrea, 26 Nov, after phase 1).

1. **Recognising chromosomes.** Chromosome names are recognised
   automatically (X, Y, XY, PAR1/PAR2, MT/M, Z, W, with or without "chr").
   Numeric codes count only when the user declares them, with
   `gt_set_chromosome_types()` (e.g. `chr_set = "human"`). Alternatives:
   - PLINK's human codes recognised in PLINK and EIGENSOFT input: rejected,
     because it is wrong for non-human data (chromosome 23 is an autosome in
     a dog, and in `lobster.bim`).
   - Declaration only: rejected, because it adds a step for every human
     dataset named with X/Y/MT.
2. **Coding.** A hemizygous genotype is stored as 0/2 (homozygous), as in
   PLINK, `.bed` files and tidypopgen's pseudohaploids, so the kernels'
   pseudohaploid path (one allele, dosage/2) carries over to phase 2. A VCF
   haploid call in a diploid individual changes from 0/1 to 0/2.
   Alternative: 0/1 with a ploidy per (individual, locus), which every
   reader and kernel would need to tell a heterozygote from a hemizygous
   ALT call.
3. **Defaults.** Analyses that pool loci and assume diploid autosomes use
   autosomal (and PAR) loci only, and say how many loci they left out.
   Alternatives: an error asking the user to filter; or using all loci with
   a warning. The escape hatch is to declare a chromosome as an autosome
   (`types = c(X = "autosome")`), so the analysis functions get no new
   argument.
4. **Scope.** The foundation and the autosome-only defaults first (phase
   1). Sex-aware statistics on X and MT come in phase 2.

## Design (phase 1)

### Chromosome types (`R/gt_chromosome_types.R`)

- **Types:** `autosome`, `X`, `Y`, `Z`, `W`, `PAR`, `MT` (MT also covers
  plastid genomes). The system (XY or ZW) follows from the types present;
  both at once is an error. The heterogametic sex is male in XY and female
  in ZW.
- **Recognised names:** matched case-insensitively, after stripping
  `chr`/`chr_`/`chromosome`/`chromosome_` (as `cast_chromosome_to_int()`
  does):
  - `X`, `Y`, `Z`, `W`;
  - `XY`, `PAR`, `PAR1`, `PAR2` as PAR;
  - `MT`, `M`, `0M` as MT.

  Numbers are never recognised.
- **`gt_set_chromosome_types(.x, types = NULL, chr_set = NULL, par = NULL)`.**
  Each call replaces the earlier declaration; a call with no declaration
  clears it.
  - `types`: a named character vector, chromosome name (as stored) to type,
    e.g. `c("23" = "X", "24" = "Y", "25" = "PAR", "26" = "MT")`.
  - `chr_set`: the number of autosomes n, as PLINK's `--chr-set`, so that
    n+1 is X, n+2 Y, n+3 PAR and n+4 MT. `"human"` is n = 22 plus
    EIGENSOFT's 90 (MT) and 91 (XY). `types` overrides it.
  - `par`: `"b36"`/`"hg18"`, `"b37"`/`"hg19"`, `"b38"`/`"hg38"` (PLINK's
    `--split-x` boundaries, for the chromosome typed X), or a data frame
    with `chromosome`, `start` and `end` (inclusive).
  - The declaration is stored as a small attribute, `chr_types`, of the
    genotypes column (a list of the types and the PAR ranges): nothing per
    locus.
- **Getters:**
  - `gt_chromosome_types(.x)` returns a tibble: chromosome, type, n_loci;
  - `loci_chromosome_types(.x)` returns one type per locus (a factor), for
    `select_loci_if()`.

  With integer chromosomes (`chromosomes_as_int`), the names come from the
  chromosome sidecar (`chromosome_sidecar_names()`).
- **Propagation:** attributes survive `filter()`, `select_loci()`,
  `gt_save()` and `gt_load()`. The functions that rebuild the column with
  `new_vctrs_bigsnp()` carry the declaration across:
  - `gt_update_backingfile()`, and so `gt_order_loci()`;
  - `rbind()`: both declarations must agree.

### Sex

- `sex` is the standard column. An internal `.gt_sex()` returns 1 (male),
  2 (female) or NA. It accepts the factor read from `.fam`,
  "male"/"female", "M"/"F" and PLINK's 1/2/0.
- The packedancestrymap reader keeps the `.ind` sex column (M/F/U).

### Coding at input

- **Haploid calls.** A haploid call in an individual of ploidy p > 1 is
  stored as 0 or p, in the cpp parser, the vcfR route, `gt_merge_vcfs()`
  and `gt_concat_vcfs()`.
- **Phased reading.** A haploid call becomes 0|0 or 1|1.
- **Ploidy guessing.** Lines on recognised X/Y/Z/W/MT chromosomes are
  skipped when guessing ploidy.

### Autosome-only defaults

An internal `.gt_autosomes_only()` returns the gen_tibble restricted to
autosomal and PAR loci, with a message, when there are loci of other types.
It looks at the chromosome levels first, so it costs O(levels) when there are
none. It is called at the start of:

- `indiv_het_obs`, `indiv_inbreeding`, `qc_report_indiv`;
- `pairwise_king`, `pairwise_ibs`, `pairwise_dist`,
  `pairwise_allele_sharing`, `pairwise_grm`;
- `gt_pca` and its `*SVD` variants;
- `gt_extract_f2`, `gt_admixture`, `gt_snmf`;
- `windows_indiv_roh`;
- the `pop_*` summaries (`pop_fst`, `pop_fis`, `pop_het_obs`,
  `pop_het_exp`, `pop_gene_div`, `pop_global_stats`), `pairwise_pop_fst`,
  `nwise_pop_pbs`, `pop_tajimas_d`, the SFS functions, and their windowed
  and jackknife variants.

`gt_pseudohaploid()` tests autosomal loci. Unchanged in phase 1: the
per-locus `loci_*` functions and `qc_report_loci`, LD clumping and pruning
(Andrea's decision of 19 Nov 2026 to treat all chromosomes alike), and
imputation.

### Sex check

`indiv_sex_check()` follows PLINK 1.9's `--check-sex`: F on non-PAR X (or Z)
loci, `max_female_f = 0.2`, `min_male_f = 0.8`, and optional Y (or W)
counts. It is tested against PLINK 1.9 output.

### Output

- **`gt_as_plink(chromosomes_as_int = TRUE)`:** X, Y, PAR and MT get n+1 to
  n+4, where n is the largest autosome number.
- **`gt_as_vcf()`:** writes haploid calls for hemizygous genotypes of
  individuals of known sex. Heterozygous haploid calls are written as
  missing.

## Phase 2 (outline)

- **Per-locus ploidy in the kernels.** Ploidy per (individual, locus class)
  goes through the per-individual lookup tables
  (`GroupedRows`/`row_class` in `src/fbm_byte_acc.h`): one table set per
  class (autosome/PAR, X, Y, MT), with loci processed class by class. The
  dataset-wide `is_pseudohaploid` flag becomes a per-class weight table.
- **Sex-aware statistics:**
  - frequencies count males once on X, as PLINK 2 `--freq`;
  - missingness leaves out female Y;
  - HWE on X follows Graffelman & Weir;
  - pi, Tajima's D, SFS and Hudson Fst on X/MT use the right sample sizes;
  - QC reports count heterozygous haploid calls, and a function erases them
    into a new backing file;
  - LD on X follows PLINK 1.9 `--ld-xchr`;
  - ROH runs on female X.

## As built (26 Nov 2026)

Differences from the design above, and details settled while building:

- **API.** Declarations take a named vector, `types = c("23" = "X")`, not
  one argument per type: arguments named `X`, `Y` or `W` would sit oddly next
  to `.x` and fail the object-name lint. `chr_set` is stored with the
  declaration, for `gt_as_plink()`.
- **Integer chromosomes.** Without declarations, only the codes whose
  names (in the chromosome sidecar) are recognised are looked up; the
  sidecar is read once with `data.table::fread()` and cached by size and
  modification time, as RAD data can have millions of codes.
- **Readers.**
  - Before the ploidies are known, a haploid ALT gets the provisional value
    254, replaced by the ploidy of its individual once it is known.
  - Haploid genotypes of individuals of ploidy > 1 are homozygotes
    everywhere, not only on sex chromosomes.
  - The ploidy guess only skips haploid genotypes on chromosomes recognised
    by name (X, Y, Z, W, MT, M, 0M), since readers know no declarations:
    `haploid_first_pop_a.vcf` (X coded as 23) still fails, as documented.
- **Sex check.** It matches PLINK 1.9 b7.11 exactly (F, SNPSEX, STATUS,
  YCOUNT) on a simulated fixture (`data-raw/validate/create_plink_sexcheck_test_data.R`).
  - From `sexcheck()` in `plink_misc.c`: F = (N − E − H)/(N − E) over
    polymorphic X loci where the individual is genotyped, E summing
    1 − 2pq, no small-sample correction.
  - From `calc_freqs_and_hwe()` in `plink_filter.c`: frequencies come from
    founders, with recorded males counting one allele (their
    heterozygotes left out), and females and unknown sex two.
  - With no polymorphic X locus, it is an error, as in PLINK.
  - Not included: the count of heterozygous haploid genotypes per
    individual (planned; left to the phase 2 QC reports).
- **Autosome-only defaults.**
  - The helper `.gt_autosomes_only()` is called at the start of every
    method (tibble, grouped and genotypes column) of the functions listed,
    so a function that calls another gives a single message.
  - A file name (`gt_admixture()`, `gt_snmf()`) passes through.
  - With no autosomal loci left, it is an error that points to the
    declaration.
- **Writers.**
  - `gt_write_text_by_locus()` takes a set of strings per individual
    (`indiv_lut`) and counts flagged strings; `gt_as_vcf()` writes runs of
    loci of the same type.
  - Individuals of unknown sex are diploid on X and haploid on Y and MT.
  - `gt_as_plink(chromosomes_as_int = TRUE)` uses n = the declared
    `chr_set`, else max(22, largest numeric autosome), PLINK's default
    being human; non-numeric autosomes get n + 4 + rank.

## Numbers

The change to the inner loop of the VCF parser (`count_alt_alleles()`)
costs nothing measurable:

- **Setup:** `gen_tibble()` on `pop1.vcf.gz`, with 100 samples and 1.02M
  biallelic SNPs (`allow_duplicates = TRUE`).
- **Where and when:** manica01, 26 Nov 2026, quiet machine (load 0.2).
- **Builds:** installed -O2 builds of `dev` (4374e70d) and of the branch
  (815611fa), A B B A A B, pinned to P-cores.

| threads | `dev` median (min–max) | branch median (min–max) |
|---|---|---|
| 1 | 2.081 s (2.080–2.088) | 2.091 s (2.090–2.107) |
| 8 | 1.314 s (1.312–1.322) | 1.330 s (1.320–1.332) |

The difference is +0.5% and +1.2%, within noise, and the backing files
are identical (md5). To reproduce: the scripts `time_parser.R` and
`run_ab.sh` in `~/project_temp/claude/tidypopgen/parser_ab/` on manica01.

## Open questions

- Phase 2: which statistics to make sex-aware first (frequencies and
  missingness on X, HWE on X, Fst on X?), and whether to add a function
  that erases heterozygous haploid genotypes (as `--set-hh-missing`).
