gt_as_plink <- function(x, file = NULL, overwrite = TRUE){

  if (is.null(file)){
    file <- bigstatsr::sub_bk(attr(x$genotypes,"bigsnp")$genotypes$backingfile,paste0(".vcf"))
  }
  if (file_ext(file)!="vcf"){
    file <- paste0(file,".vcf")
  }
  if (!overwrite && file.exists(file)){
    stop("file ", file," already exists; use 'overwrite=TRUE' to overwrite it")
  }

  # if chunk is null, get the best guess of an efficient approach
  if (is.null(chunk_size)){
    chunk_size <- bigstatsr::block_size(nrow(x))
  }

  # set up chunks
  chunks_vec <- c(
    rep(chunk_size, floor(count_loci(x) / chunk_size)),
    count_loci(x) %% chunk_size
  )
  chunks_vec_index <- c(1,cumsum(chunks_vec))

  # generate the header
  vcf_header <- c("##fileformat=VCFv4.3",
                  paste0("##fileDate=",format(Sys.time(), "%Y%m%e")),
                  paste0("##source=tidypopgen_v",packageVersion("tidypopgen")))
  # create copy of loci table
  loci_tbl <- show_loci(x)
  # reorder chromosomes levels in the order in which they appear
  loci_tbl$chromosome <- forcats::fct_inorder(loci_tbl$chromosome)
  # get max position for each chromosome
  chromosome_summary <- loci_tbl %>% group_by(chromosome) %>% summarise(max = max(.data$position))
  for (chrom_i in 1:nrow(chromosome_summary)){
    vcf_header <- c(vcf_header,
                    paste0("##contig=<ID=", chromosome_summary$chromosome[chrom_i],
                           ",length=",chromosome_summary$max[chrom_i]+1,
                           ">"))
  }
  vcf_header <- c(vcf_header, '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">')
  vcf_header <- c(vcf_header, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  vcf_header <- c(vcf_header, paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                                     paste(x$id, collapse="\t")))
  # write it to file
  write.table(vcf_header, file = file, quote = FALSE, row.names = FALSE,
              col.names = F)
  # iterate over genotypes in chunks and append to the vcf body
  for (chunk_i in seq_along(chunks_vec)) {
    genotypes_matrix <- t(show_genotypes(x,
                                         loci_indices =
                                           chunks_vec_index[chunk_i]:chunks_vec_index[chunk_i+1]))
    genotypes_matrix[genotypes_matrix==0] <- "0/0"
    genotypes_matrix[genotypes_matrix==1] <- "0/1"
    genotypes_matrix[genotypes_matrix==2] <- "1/1"
    genotypes_matrix[is.na(genotypes_matrix)] <- "./."
    # subset loci to this chunk
    loci_sub <- show_loci(x)[chunks_vec_index[chunk_i]:chunks_vec_index[chunk_i+1],]
    # add the other columns needed for the
    loci_cols <- c("chromosome", "position", "name", "allele_ref", "allele_alt")
    loci_sub <- loci_sub %>% select(any_of(loci_cols)) %>%
      mutate(qual=".",filter=".", info="PR", format = "GT")
    loci_sub[is.na(loci_sub)] <- "."
    genotypes_matrix <- cbind(loci_sub, genotypes_matrix)
    # append table to previous chunk
    write.table(genotypes_matrix, file = file, quote = FALSE, append = TRUE,
                col.names = FALSE, row.names = FALSE)
  }

}

# TODO for tests
# read pop_a_vcf
# rewrite it as a vcf
# re-read it and check that we got back what we started with


vcf_path <- system.file("extdata/pop_a.vcf", package = "tidypopgen")
pop_a_vcf_gt <- gen_tibble(vcf_path, quiet=TRUE,backingfile = tempfile())

x <- pop_a_vcf_gt
chunk_size <- 10
