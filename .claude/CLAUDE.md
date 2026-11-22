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
