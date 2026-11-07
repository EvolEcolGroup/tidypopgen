# Managing the manual review of `codereview`

This note describes how to check the commits on `codereview` (listed in
[code_review_2026-10.md](code_review_2026-10.md)) before they are merged
into `dev`. It is written for collaborators and for coding agents: an agent
given this file should follow it literally, and stop wherever it says
**ask the maintainer**.

## Rules

1. **The report is the record.** Every commit has a checkbox in the report.
   A box is ticked only once its commit has passed the checks below, and
   review status is not tracked anywhere else.
2. **One commit at a time.** Start with
   [Changes in behaviour](code_review_2026-10.md#changes-in-behaviour-to-check-before-merging),
   then go through [Fixes by area](code_review_2026-10.md#fixes-by-area).
   Claim a section by adding a line `Reviewer: <your initials>` under its
   heading (not in the heading, which would break the links to it), so
   that two people do not check the same commits.
3. **Never rewrite history.** No rebase, amend, squash or force push on
   `codereview`: the report refers to commits by their hash. Every change,
   including a revert, is a new commit on top.
4. **Speed matters as much as correctness.** A fix that makes a hot path
   slower is not accepted as it is (see [Performance](#performance)).

## Checking a commit

For a commit `<hash>`:

1. **Read it** with `git show <hash>`: the message, the code, the test and
   the `NEWS.md` entry. Check that the problem is real, that the fix matches
   the documentation, and that it is the smallest sensible change.
2. **Check that the test catches the bug.** In a separate worktree, put back
   the old code but keep the new tests. The test must fail, and pass again
   on the commit itself:

   ```sh
   git worktree add ../tpg_check <hash>
   cd ../tpg_check
   # old version of the files the commit modified, without the files it added
   git diff --name-only --diff-filter=M <hash>~1 <hash> -- . ':!tests' | xargs -r git checkout <hash>~1 --
   git diff --name-only --diff-filter=A <hash>~1 <hash> -- . ':!tests' | xargs -r rm
   Rscript -e 'devtools::test(filter = "<test file name without test_>")'  # must fail
   git checkout <hash> -- .
   Rscript -e 'devtools::test(filter = "<test file name without test_>")'  # must pass
   ```

   Skip this step for commits listed under "No new test" in the report's
   Verification section.
3. **Behaviour changes** (commits in "Changes in behaviour"): whether the new
   behaviour is acceptable for users is a design decision. **Ask the
   maintainer.**
4. **Performance:** if the commit touches a hot path, benchmark it (next
   section).
5. **Record the outcome** in the report (see
   [Recording the outcome](#recording-the-outcome)).

## Performance

Hot paths are the C++ code in `src/` (the VCF and packedancestry readers,
and the per-locus and pairwise kernels), the file readers behind
`gen_tibble()`, and any R code that loops over loci or individuals or reads
the backing file.

For these, compare the commit with `dev` before ticking its box:

- Install both versions in separate libraries (see the top of
  `data-raw/benchmarking/benchmark_vcf_parser.R` for the commands), and
  run each version in its own R session, at least three times each, in a
  balanced order (e.g. A B B A A B): a version that always runs in the
  same position can look faster or slower than it is. Keep the machine
  from sleeping (e.g. with `systemd-inhibit`), as elapsed times include
  the time spent suspended.
- Use data at the scale users have, not the test data. For VCFs, that means
  50,000 variants × 1,000 samples, the same with one sample that is `.` on
  every variant, 10 million variants × 20 samples (plain and gzipped), and
  tetraploid and mixed-ploidy files. `benchmark_vcf_parser.R` creates these
  files and times the cpp parser.
- The outputs must be identical (e.g. the md5 of the `.bk` file), unless
  the commit is meant to change them.
- **Accept** only if there is no slowdown beyond run-to-run noise (about
  5%). Otherwise **ask the maintainer**, with the timings. The usual remedy
  is a new commit that keeps the fix and restores the speed, as `7081979`
  did for the VCF fixes.

## Recording the outcome

Edit the commit's entry in the report:

- **Accepted:** tick the box (`[x]`) and add your initials at the end of
  the entry, e.g. `— AM`.
- **Accepted after a follow-up fix:** add the follow-up commit as a new entry
  in the same section (same format: hash · functions — problem *[impact]*),
  mention it in the original entry ("reimplemented in `<new hash>`"), and
  tick both once the follow-up has been checked.
- **Rejected:** revert it with `git revert <hash>`. Keep its entry, add
  "Reverted in `<revert hash>`: <reason>", and add the revert as an entry.
- **Decision needed:** leave the box empty, add "**Decision needed:**
  <question>", and **ask the maintainer**.

Commit changes to the report separately from code changes, so that every
code commit can be checked on its own. After any code change, re-run the
checks in the report's Verification section and update it.

## Instructions for agents

- Work only on the commits or section you were given. For each commit,
  report its hash, the outcome (accepted / problem found / decision needed)
  and the evidence (test output, benchmark timings).
- You may run tests and benchmarks, write follow-up commits, and edit the
  report entries of the commits you checked.
- Do not tick a box for a check you could not run (e.g. a dependency is
  missing): say so in the entry instead.
- Do not push, open or merge pull requests, or merge into `dev` without the
  maintainer's explicit go-ahead.
- Stop and **ask the maintainer** when a test does not fail on the old
  code, when a benchmark is slower, when outputs differ from `dev` where no
  change was intended, or when a fix would change the API.
- Use worktrees and separate libraries, so that the maintainer's checkout
  is never modified, and delete large generated files when done.

## Before merging into `dev`

1. Every box is ticked, or has a recorded decision.
2. The checks in the report's Verification section have been re-run on the
   tip of `codereview`, with the package installed, and the table updated.
3. The [Findings not fixed](code_review_2026-10.md#findings-not-fixed) have
   been decided or turned into issues, and the optional branch
   `codereview-bigsnp-rds` has been merged or dropped.
4. A pull request from `codereview` to `dev` links the report. Copying a
   section of the report into the pull request makes its boxes clickable.
