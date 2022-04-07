#' Get Codons From Gene Sequence
#'
#' Group a gene sequence into codons (triplets) that can be used for subsequent translation. If any elements of the character vector have length >1, they are trimmed to the first base. Primarily used by id.snps and id.indels functions.
#' @param gene.seq Character vector containing gene sequence
#' @keywords codons
#' @export
#' @examples
#' get_codons()

get_codons <- function(gene.seq) {
  CodonStartPositions <- 0:(length(gene.seq)/3-1) * 3 + 1
  CodonEndPositions <- 1:(length(gene.seq)/3) * 3
  #this function will be applied to the ALT column when run from call_mutations()
  #entries in this column will have length > 1 in the case of insertions
  #insertions will be dealt with separately, so here, trim those to the first nucleotide
  first_base <- stringr::str_sub(gene.seq, start = 1L, end = 1L)
  gene_seq <- first_base %>% unlist %>% paste0(collapse = "")
  codons <- stringr::str_sub(string = gene_seq, start = CodonStartPositions, end = CodonEndPositions)
  rep(codons, each = 3)
}
