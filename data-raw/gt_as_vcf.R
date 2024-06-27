# if chunk is null, get the best guess of an efficient approach
if (is.null(chunk_size)){
  chunk_size <- bigstatsr::block_size(nrow(x))
}

# set up chunks
chunks_vec <- c(
  rep(chunk_size, floor(count_loci(x) / chunk_size)),
  count_loci(x) %% chunk_size
)
chunks_vec_index <- c(1, chunks_vec)

# generate the header
vcf_header <- c("##fileformat=VCFv4.3",
                paste0("##fileDate=",format(Sys.time(), "%Y%m%e")),
                paste0("##source=tidypopgen",packageVersion("tidypopgen")))
chromosome_summary <- show_loci(x) %>% group_by(chromosome) %>% summarise(max = max(.data$position))
for (chrom_i in 1:nrow(chromosome_summary)){
  vcf_header <- c(vcf_header,
                  paste0("##contig=<ID=", chromosome_summary$chromosome[chrom_i],
                         ",length=",chromosome_summary$max[chrom_i]+1,
                  ">"))
}
vcf_header <- c(vcf_header, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
vcf_header <- c(vcf_header, paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t",
                                   paste(x$id, collapse="\t")))
# write it to file

# iterate over genotypes in chunks and append to the vcf body
for (chunk_i in seq_along(chunks_vec)) {
  genotypes_matrix <- show_genotypes(x,loci_indices = chunks_vec_index[chunk_i]:chunks_vec_index[chunk_i]+1)
}


  ##fileformat=VCFv4.3
  ##fileDate=20090805
  ##source=myImputationProgramV3.1
  ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
  ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
  ##phasing=partial
  ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
