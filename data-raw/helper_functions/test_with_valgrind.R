R -d valgrind --vanilla -e "devtools::load_all(); testthat::test_file('tests/testthat/test_pairwise_pop_fst.R')"


R -d "valgrind --tool=memcheck --leak-check=full" --vanilla -e "devtools::load_all(); testthat::test_file('tests/testthat/test_augment_loci')"