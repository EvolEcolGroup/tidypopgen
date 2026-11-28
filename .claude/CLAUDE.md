# Notes for Claude on tidypopgen

## Never benchmark -O0 builds of the C++ code

- `pkgload::load_all()` / `devtools::load_all()` compile `src/` with pkgbuild's
  debug flags (`-UNDEBUG -Wall -pedantic -g -O0`). The kernels then run 2-3x
  slower than an installed (-O2) build, so never time or profile compiled code
  from a plain `load_all()`.
- To time a checkout: `pkgbuild::clean_dll(path)`, then
  `pkgbuild::compile_dll(path, debug = FALSE)`, then
  `pkgload::load_all(path, compile = FALSE)`. After editing anything in `src/`,
  run `compile_dll(path, debug = FALSE)` again before loading. Otherwise
  `load_all()` recompiles the changed files at -O0 and `src/tidypopgen.so`
  silently mixes -O0 and -O2 objects. `compile_dll(force = TRUE)` does not
  rebuild up-to-date `.o` files: use `clean_dll()` first.
- For a baseline, install the reference commit into its own library:
  `git archive <ref> | tar -x -C <dir>`, then
  `R CMD INSTALL --library=<lib> <dir>`, and load it in a separate R process
  with `library(tidypopgen, lib.loc = <lib>)`. Never `load_all()` a checkout
  that other processes or subagents may be compiling at the same time.
- To check a build: `objdump -d src/tidypopgen.so | grep -c 'call.*vectorI.*ixEm'`
  greater than 0 means -O0 code is present. Per object,
  `readelf --debug-dump=info src/<file>.o | grep -m1 DW_AT_producer` shows the
  flags (the last `-O` wins).
- When briefing subagents that build or benchmark this package, pass these
  rules on explicitly.

## `[skip ci]` skips CI for the whole push, and for squash merges

- GitHub skips every workflow of a push (and of a PR) when the head commit's
  message contains `[skip ci]` anywhere. A code commit pushed together with
  a later `[skip ci]` notes commit is never checked: the TGENO reader
  (`8284e14d`) reached `dev` that way, with a test that failed on Windows
  (fixed in #248).
- A squash merge lists the messages of all the squashed commits, so a PR
  with one `[skip ci]` commit gives a squash commit that skips CI on `dev`
  (#246), and the `dev` into `main` PR of a release skips pkgdown on `main`,
  so the site is not rebuilt: start it by hand
  (`gh workflow run pkgdown.yaml --ref main`), or remove `[skip ci]` from
  the squash message.
- So: push code commits on their own (or end the push with a commit without
  `[skip ci]`), and check that CI ran (`gh run list --commit <sha>`).
- `main` takes PRs only through a merge queue (squash; no required checks):
  `gh pr merge` fails there ("merge strategy set by the merge queue", and
  auto-merge is off). Enqueue with
  `gh api graphql -f query='mutation($id:ID!){enqueuePullRequest(input:{pullRequestId:$id}){mergeQueueEntry{position state}}}' -f id=$(gh pr view <n> --json id -q .id)`.
  Its squash message lists all the commits, so a release from `dev` skips
  CI and pkgdown on `main`. Then merge `main` back into `dev` (no file
  changes), so that `main` stays an ancestor of `dev`.
