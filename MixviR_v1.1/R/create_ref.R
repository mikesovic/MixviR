#' Create MixVir reference genome object
#'
#' Uses a fasta genome (single chromosome) and bed file defining features of interest (genes) to create a data frame used as a reference to call variants/mutations from a sample.
#' @param genome fasta formatted genome file (assumes only one chromosome)
#' @param feature.bed bed file (4 columns, tab delimited, no column names: chr, start, end, feature_name) with info on features of interest (open reading frames to translate)
#' @keywords reference
#' @export
#' @examples
#' create_ref()

create_ref <- function(genome, feature.bed) {

  features <- readr::read_tsv(feature.bed, col_names = FALSE)
  names(features) <- c("chrm", "start", "end", "feature")

  features <- features %>%
    dplyr::group_by(r=dplyr::row_number()) %>%
    dplyr::mutate("pos" = list(seq.int(from = start, to = end))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-r) %>%
    tidyr::unnest(cols = c(pos))


  sequence <- Biostrings::readDNAStringSet(genome)

  ref_genome <- data.frame("chr" = character(),
                           "pos" = integer(),
                           "ref_base" = character()
                           )

  for (i in length(sequence)) {
    chr_seq <- as.character(sequence[[i]]) %>%
      stringr::str_split(pattern = "") %>%
      unlist()
    chr_seq <- data.frame("chr" = names(sequence)[i],
                          "pos" = 1:length(sequence[[i]]),
                          "ref_base" = chr_seq)
    ref_genome <- dplyr::bind_rows(ref_genome, chr_seq)
  }

  merged <- dplyr::left_join(x = ref_genome,
                             y = features,
                             by = "pos") %>%
    dplyr::select(-chrm)

  feature_positions <- merged %>% dplyr::filter(!is.na(feature))

  split_to_codons <- function(x) {
    StartVec <- 0:(length(x)/3-1) * 3 + 1
    EndVec <- 1:(length(x)/3) * 3
    gene_seq <- x %>% unlist %>% paste0(collapse = "")
    codons <- stringr::str_sub(string = gene_seq, start = StartVec, end = EndVec)
    rep(codons, each = 3)
  }

  feature_positions <- feature_positions %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate("ref_codon" = split_to_codons(ref_base))

  codons <- Biostrings::DNAStringSet(feature_positions$ref_codon)
  aas <- Biostrings::translate(codons, no.init.codon = TRUE)
  aas <- as.character(aas, use.names = FALSE)

  feature_positions <- feature_positions %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate("gene_aa_position" = rep(1:(dplyr::n()/3), each = 3)) %>%
    dplyr::mutate("ref_AA" = aas)

  nonfeature_positions <- merged %>%
    dplyr::filter(is.na(feature)) %>%
    dplyr::mutate("ref_codon" = NA) %>%
    dplyr::mutate("ref_AA" = NA) %>%
    dplyr::mutate("gene_aa_position" = NA)

  all_ref <- dplyr::bind_rows(nonfeature_positions,
                       feature_positions) %>%
    dplyr::arrange(genomic_pos)

  all_ref$feature[which(is.na(all_ref$feature))] <- "non-genic"

  all_ref <- all_ref %>% tidyr::unite(col = "ref1",
                                      feature, ref_AA,
                                      sep = "_",
                                      remove = FALSE) %>%
    tidyr::unite(col = "ref_identity",
          ref1, gene_aa_position,
          sep = "",
          remove = FALSE) %>%
    dplyr::select(-ref1) %>%
    dplyr::rename("gene" = "feature")

  all_ref <- all_ref %>%
    dplyr::select(genomic_pos,
                 ref_base,
                 gene,
                 ref_codon,
                 ref_AA,
                 gene_aa_position,
                 ref_identity
    )

  all_ref
}
