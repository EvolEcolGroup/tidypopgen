# CRAN allows at most 2 cores: limit the cores used by all tests. With
# option mc.cores set to 2, bigstatsr::nb_cores() (the default n_cores of most
# functions) returns 1; tests that set n_cores themselves use at most 2.
withr::local_options(mc.cores = 2, .local_envir = teardown_env())
data.table::setDTthreads(2)
if (rlang::is_installed("RhpcBLASctl")) {
  RhpcBLASctl::blas_set_num_threads(2)
  RhpcBLASctl::omp_set_num_threads(2)
}
