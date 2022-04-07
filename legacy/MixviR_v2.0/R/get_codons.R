#' Get Codons From Gene Sequence
#'
#' Group a gene sequence into codons (triplets) that can be used for subsequent translation. If any elements of the character vector have length >1 (insertions in ALT column of VCF), they are trimmed to the first base. Used by id.snps, id.indels, and create_ref functions.
#' @param gene.seq Character vector containing gene sequence, with length of vector equal to length of sequence.
#' @keywords codons
#' @examples
#' get_codons()

get_codons <- function(gene.seq) {
  CodonStartPositions <- 0:(length(gene.seq)/3-1) * 3 + 1
  CodonEndPositions <- 1:(length(gene.seq)/3) * 3
  first_base <- stringr::str_sub(gene.seq, start = 1L, end = 1L)
  gene_seq <- first_base %>% 
    unlist() %>% 
    paste0(collapse = "")
  codons <- stringr::str_sub(string = gene_seq, start = CodonStartPositions, end = CodonEndPositions)
  rep(codons, each = 3)
}
