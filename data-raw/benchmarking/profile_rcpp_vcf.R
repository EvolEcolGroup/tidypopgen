library(Rgperftools)
vcf_path <- "../../Downloads/vcf_test/pileup_Jar_burial.pileUp_Choskar.fixREF.vcf.gz"
profile_file <- file.path(tempdir(),"vcf_profile3.out")
start_profiler(profile_file)
pop_a_vcf_fast_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(), parser="cpp")
stop_profiler()


## from bash
~/go/bin/pprof --http=localhost:8080 ~/git/tidypopgen/src/tidypopgen.so /tmp/RtmpTmt48o/vcf_profile2.out




##############
# profile in R
#
Rprof()
pop_a_vcf_fast_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(), parser="cpp")
Rprof(NULL)
summaryRprof()


##
library(tictoc)
tic()
pop_a_vcf_fast_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile(), parser="cpp")
toc()
