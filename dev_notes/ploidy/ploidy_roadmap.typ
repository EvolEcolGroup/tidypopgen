#import "@local/report:0.1.0": *
#show: report.with(
  title: "Arbitrary ploidy in tidypopgen",
  subtitle: "Roadmap, from the state of v1.0.0",
  abstract: [
    tidypopgen stores genotypes as allele dosages with a ploidy for each
    individual, so the data model already copes with haploids, polyploids and
    mixed-ploidy samples. The analyses do not: of 59 calls tried on tetraploid
    data at v1.0.0, 34 stop with an error and 24 run, three of them warning
    that they are untested and two returning numbers that assume diploidy. This note records where
    each function stands (from an audit script that calls every function on
    simulated data), sets out the order in which to close the gaps, and gives
    the estimator and the validation planned for each step.
  ],
  toc: true,
)

= Where we are <sec-state>

A `gen_tibble` keeps one byte per genotype, the dosage of the alternate allele,
in a file-backed matrix. The genotype column carries two attributes: `ploidy`,
a single label (the ploidy if all individuals share it, `0` for mixed ploidy
and `-2` when some individuals are pseudohaploid), and `fbm_ploidy`, the ploidy
of every individual in the backing file. `gen_tibble()` reads polyploid and
mixed-ploidy VCF files and matrices, and `indiv_ploidy()` and `show_ploidy()`
expose the two attributes. Nothing in this model is specific to diploids, and
the fast byte accessor in `src/fbm_byte_acc.h` already takes a lookup table for
each individual, which is what a ploidy-aware kernel needs.

@tab-audit shows what happens when each analysis is called on simulated data
with three populations of twelve individuals and 600 loci, for five ploidy
configurations. It is produced by `audit.R`, next to this file, and is the
measure of progress for the whole roadmap: the aim is a table with no
unexplained `error`, and no `ok` that hides a diploid assumption.

#figure(
  tsv-table(read("audit.tsv")),
  caption: [Outcome of each call on simulated data, by ploidy (`audit.R`).
    `ok` means only that the call returned without a condition, not that the
    result is right for that ploidy; see the text. `no_pkg` marks a call that
    needs a package that was not installed where the audit ran (admixtools).
    This is the table for the working branch (see Progress); the one for
    v1.0.0, which the rest of this section describes, is kept in
    `audit_v1.0.0.tsv`, and the messages in `audit_long.tsv`.],
) <tab-audit>

The functions fall in three groups. The first is ploidy-aware by design and
tested: reading data, missingness, ungrouped allele frequencies, `loci_pi()`,
`indiv_het_obs()`, `pop_het_obs()`, Tajima's D and the site frequency spectra
with their windowed and jackknifed versions, and `gt_as_genind()`. The second
stops with "this function only works on diploid data": expected
heterozygosity, $F_"IS"$, $F_"ST"$, the global statistics, Hardy-Weinberg
tests, identity by state, allele sharing, the GRM, KING, inbreeding, LD
clumping, the locus QC report, `rbind()` and most exporters. The third group
is small, but it is the one to worry about: imputation, the PCA functions and
what is built on them, and runs of homozygosity have no ploidy check at all.
`gt_pca_autoSVD()` and `windows_indiv_roh()` return a result for every ploidy,
computed as if the dosages were diploid, and all the PCA functions do the same
for haploids. `gt_pca_partialSVD()` and `gt_pca_randomSVD()` stop on polyploids
only by accident, because the binomial scaling of bigsnpr has no value for an
allele frequency above one, and imputation stops with a message about byte
codes.

#caution[
  An `ok` in the polyploid columns for a function in the third group is a bug,
  not a feature. The first piece of work is to make every function either
  correct or explicit, before any estimator is generalised.
]

= Principles

Three rules keep the work honest. First, the unit of observation is the allele
copy: sample sizes are sums of ploidy over individuals with a genotype, never
twice the number of individuals. Every statistic that depends on the data only
through allele counts (frequencies, $pi$, Hudson's $F_"ST"$, $f_2$, the SFS,
Tajima's D) then generalises without new theory. Second, a general kernel must
return exactly the diploid result when every ploidy is two, and the existing
diploid kernels stay as fast paths behind the same R function; no user-facing
function is split by ploidy. Third, where a statistic has no agreed definition
beyond diploids (KING, PLINK's runs of homozygosity, exact Hardy-Weinberg tests
under polysomic inheritance), the function says so and points to the
alternative, rather than inventing an estimator.

Validation has three layers, in increasing cost. Pure R reference
implementations on the dosage matrix, as already written for
`test_loci_pi_polyploid.R`. An exact invariant for allele-count statistics:
merging pairs of diploid individuals of the same population into tetraploids
(adding their dosages) leaves the allele counts unchanged, so the tetraploid
result must equal the diploid one to machine precision. And comparison with
the polyploid packages, with fixtures generated once in `data-raw/validate` as
for easySFS and pixy: adegenet and hierfstat through `gt_as_genind()`, StAMPP
(mixed-ploidy $F_"ST"$, Nei's distance and GRM), AGHmatrix (VanRaden
relationship matrices for autopolyploids) and GenoDive ($rho$ and polyploid
AMOVA).

= The roadmap

== Phase 0: foundations

One internal description of ploidy replaces the four overlapping helpers
(`stopifnot_diploid()`, `stopifnot_dip_pseudo()`, `is_diploid_only()`,
`stopifnot_no_pseudohaploid()`). A function `ploidy_class()` returns one of
`haploid`, `diploid`, `pseudohaploid`, `polyploid` or `mixed`, computed from
the individuals currently in the tibble, so the label can no longer go stale
after subsetting (`.gt_refresh_ploidy()` exists only to patch that). A single
guard, `check_ploidy(.x, supported = ...)`, produces the same error everywhere,
naming the calling function, the ploidy found and what to use instead. Every
exported function that touches genotypes calls it, which closes the third
group of @sec-state at once: PCA, imputation and runs of homozygosity become
explicit errors until their phase arrives. `qc_report_indiv()`, which runs on
every ploidy, is read through for diploid assumptions at the same time.

The test infrastructure comes with it: a helper that simulates polyploid and
mixed-ploidy `gen_tibble`s with population structure and missing data, the
merging invariant above as a reusable expectation, and `audit.R` turned into a
snapshot test, so that a change in the ploidy support of any function shows up
in review. The generated table also feeds the documentation (phase 6).

== Phase 1: allele frequencies and gene diversity

Grouped `loci_alt_freq()` and `loci_maf()` currently stop for polyploids and
tell the user to loop over groups; the grouped kernel is already general apart
from the pseudohaploid multiplier, so this is a small change, and it lets the
ungrouped polyploid path drop its "needs a proper unit test" warning.
`pop_het_exp()` and `pop_gene_div()` follow. Writing Nei's correction for $M$
allele copies, $M\/(M-1)$ in place of $2n\/(2n-1)$, would assume random mating
within individuals, and would not reduce to the diploid estimator, which
corrects for observed heterozygosity. The estimator used instead is the
proportion of pairs of allele copies from two different individuals that
differ: with $A = sum_i d_i$ and $M = sum_i k_i$ over the individuals with a
genotype,
$ hat(h)_s = 2 (A (M - A) - sum_i d_i (k_i - d_i)) / (M^2 - sum_i k_i^2). $
It never compares copies within an individual, so it is unbiased whatever the
inbreeding and the mix of ploidies, it is Nei's (1987) estimator for diploids
and $n\/(n-1) (1 - sum p^2)$ for haploids. The same four sums are what Nei's
$F_"ST"$ and $F_"IS"$ need in phase 2. `qc_report_loci()`
reports MAF and missingness for any ploidy and leaves out the Hardy-Weinberg
test, as it already does for pseudohaploids.

== Phase 2: differentiation between populations

Hudson's $F_"ST"$ (Bhatia et al. 2013) uses only allele counts and sample
sizes in allele copies, and the kernel already handles pseudohaploids, so it
is extended first; `windows_pairwise_pop_fst()`, `jackknife_pairwise_pop_fst()`
and the PBS functions inherit it. Nei's 1987 estimator needs observed
heterozygosity, for which the gametic heterozygosity already used by
`pop_het_obs()` on polyploids is the natural replacement. Weir and Cockerham's
1984 estimator is tied to diploid genotypes through its within-individual
variance component; for polyploids the literature (Ronfort et al. 1998;
Meirmans and Liu 2018) recommends $rho$, which is independent of ploidy and of
the rate of double reduction, so `method = "rho"` is added and `WC84` stays
diploid-only with a message saying so.

The Weir and Goudet (2017) family (`pairwise_allele_sharing()`, `pop_fst()`,
`pop_fis(method = "WG17")`, `indiv_inbreeding()` and `pairwise_grm()`) rests on
one quantity, the probability that two alleles drawn from two individuals are
identical in state. With dosages $d$ and ploidies $k$ this is
$(d_i d_j + (k_i - d_i)(k_j - d_j)) \/ (k_i k_j)$ between individuals, and the
same expression for draws without replacement within an individual. Once the
allele sharing matrix is general, the functions built on it need no further
change. `pop_global_stats()` and `pop_fis(method = "Nei87")` close the phase.
`gt_extract_f2()` is checked against ADMIXTOOLS 2, which already takes a
ploidy for each individual.

== Phase 3: relationships between individuals

`pairwise_ibs()` is redefined on the frequency scale, one minus the mean of
$|d_i\/k_i - d_j\/k_j|$, which is the PLINK proportion for diploids; the count
outputs remain diploid-only. `pairwise_grm()` gains VanRaden's estimator with
variance $k p (1 - p)$ for autopolyploids, validated against AGHmatrix. KING's
robust estimator counts diploid genotype classes and stays diploid-only;
`qc_report_indiv()` and `filter_high_relatedness()` use the allele sharing
kinship for other ploidies. Method-of-moments relatedness for polyploids
(Huang et al. 2014) is a possible later addition.

== Phase 4: multivariate analyses and imputation

For a single ploidy $k$, PCA needs only a scaling function with centre $k p$
and scale $sqrt(k p (1 - p))$, passed to the bigstatsr SVD in place of
`snp_scaleBinom()`. Mixed ploidy needs dosages on the frequency scale with a
weight for each individual, which bigstatsr's column scaling cannot express:
the choice is between writing a scaled copy of the matrix and matrix-free
products in C++ for a partial SVD, to be decided on a benchmark. Imputation
needs its own byte codes, because bigsnpr's imputed codes assume dosages of 0,
1 and 2; mode, mean and random imputation are simple to write over the byte
accessor. DAPC, clustering and pcadapt follow the PCA. `loci_ld_clump()`
works on correlations of dosages and should extend to any single ploidy once
the imputed codes exist. Runs of homozygosity are a diploid concept as
implemented, and stay so. Of the external programs, ADMIXTURE is
diploid-only, and sNMF takes a ploidy, which `gt_snmf()` should take from the
data.

== Phase 5: input, output and combining data

`gt_as_vcf()` writes polyploid genotypes (the reader already handles them),
`rbind()` combines tibbles of any ploidy, and `gt_as_genlight()` and
`gt_from_genlight()` use the ploidy slot of genlight. PLINK, hierfstat and the
LEA geno format are diploid formats and keep their errors. A STRUCTURE export
would reach most of the polyploid clustering programs; polysat and StAMPP are
reachable through genind and genlight.

== Real data

All the tests so far use simulated dosages. A real polyploid dataset is needed
for what simulations miss (dosage errors, missingness that differs between
ploidies, the conventions of polyploid genotype callers), and for the vignette
and the benchmarks. The ideal one is public, with a licence that allows a
subset in the package, called with a polyploid-aware caller, with several
populations, both diploids and polyploids of one species, and published
statistics to compare with. The diploid and autotetraploid populations of
_Arabidopsis arenosa_ are the first candidate to check; the list is in
`README.md`. Choosing it does not depend on any phase, and is best done early.

== Phase 6: documentation

A vignette on polyploid and mixed-ploidy data, on that dataset; a "Ploidy"
section in the help of every function, generated from the same table as
@tab-audit; and the support table as a pkgdown article.

== Later: ploidy that varies along the genome, and uncertain dosages

Two extensions change the data model and are kept apart. Sex chromosomes,
haplodiploid species and organelles need a ploidy for each individual and
locus; the natural design is a ploidy class for each chromosome combined with
a sex for each individual, which the per-individual lookup tables can follow by
switching tables at chromosome boundaries. Polyploid genotype calls are often
uncertain, and updog and polyRAD output posterior mean dosages or genotype
likelihoods; supporting them means storing something other than an integer
dosage, and is out of scope until the integer case is complete.

= Progress

The work is cut into small branches off `dev`, one for each reviewable step:
`ploidy_guards` (the classifier and the guard, on which the others are
stacked), `ploidy_guard_unguarded` (PCA, imputation, pcadapt and runs of
homozygosity), `ploidy_loci_freq`, `ploidy_qc_report_loci` (on top of
`ploidy_loci_freq`), `ploidy_pop_het_exp` and `ploidy_hudson_fst`; this note is
on `ploidy_roadmap`. @tab-audit is the table with all of them applied to
v1.0.0: on tetraploid data, 31 calls run and 27 stop, against 24 and 34 at
v1.0.0, and none of those that run assumes diploidy.

#solid Phase 0 is in place apart from the snapshot test. `.ploidy_class()` and
`.check_ploidy()` are in `R/gt_helper_functions.R`, the two older guards call
them, and the PCA functions, imputation, pcadapt and runs of homozygosity
check ploidy. The test helpers are in `tests/testthat/helper_ploidy.R`:
simulation for any ploidy, reference counts and gene diversity, and the merging
invariant.

#solid Phase 1 is done for a single population or several: grouped
`loci_alt_freq()` and `loci_maf()`, `pop_het_exp()` and `pop_gene_div()` (new
kernel `src/grouped_gene_div.cpp`), and `qc_report_loci()`. The new tests found
that `qc_report_loci()` halved the allele frequencies of true haploids, which
is fixed.

#solid Phase 2 has started with Hudson's $F_"ST"$, which
`pairwise_pop_fst()`, `windows_pairwise_pop_fst()`, the jackknife and the PBS
functions now compute for any ploidy. It is tested against the formula of
Bhatia et al. applied to allele counts, and with the merging invariant.
Polyploid data go through the kernel's slower path that sums doubles rather
than lookup tables, which is the place to look if it proves slow on whole
genomes.

#caveat The global estimate of `pop_het_exp(include_global = TRUE)` is still
diploid-only; it belongs with `pop_global_stats()` in phase 2. The validation
so far is against R reference implementations and the diploid kernel, not yet
against other packages.

= Open decisions

The inheritance mode (polysomic for autopolyploids, disomic for allopolyploids)
is not recorded in a `gen_tibble`. It matters only for expectations about
genotype frequencies, so for Hardy-Weinberg tests and $F_"IS"$; the proposal is
to treat dosages as unordered allele copies throughout, and to leave
`loci_hwe()` diploid-only. Whether `pairwise_pop_fst()` should default to
Hudson for polyploids (it is the only method valid for every ploidy) is to be
agreed with the other maintainers, as is the question of whether this work
goes upstream as one pull request for each phase.
