#' Create MixVir reference genome object
#'
#' Uses a fasta genome (single chromosome) and bed file defining features of interest (genes) to create a data frame used as a reference to call variants/mutations from a sample.
#' @param genome fasta formatted genome file (assumes only one chromosome)
#' @param feature.bed bed file (4 columns, tab delimited, no column names: chr, start, end, feature_name) defining features of interest (open reading frames to translate)
#' @keywords reference
#' @export
#' @examples
#' create_ref()

create_ref <- function(genome, feature.bed) {

  features <- readr::read_tsv(feature.bed, col_names = FALSE)
  names(features) <- c("chrm", "start", "end", "GENE")

  #create df that includes row for every position in each feature/gene
  features <- features %>%
    dplyr::group_by(r=dplyr::row_number()) %>%
    dplyr::mutate("POS" = list(seq.int(from = start, to = end))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-r) %>%
    tidyr::unnest(cols = c(POS))

  #read in the genome
  sequence <- Biostrings::readDNAStringSet(genome)

  ref_genome <- data.frame("CHR" = character(),
                           "POS" = integer(),
                           "REF_BASE" = character()
                           )

  for (i in length(sequence)) { #loop over each chromosome in the genome
    chr_seq <- as.character(sequence[[i]]) %>%
      stringr::str_split(pattern = "") %>%
      unlist()
    chr_seq <- data.frame("CHR" = names(sequence)[i],
                          "POS" = 1:length(sequence[[i]]),
                          "REF_BASE" = chr_seq)
    ref_genome <- dplyr::bind_rows(ref_genome, chr_seq)
  }

  #to analyzed multiple chromosomes, need to merge on CHR+POS
  #right now, only merging on POS, so assumes one chr only
  merged <- dplyr::left_join(x = ref_genome,
                             y = features,
                             by = "POS") %>%
    dplyr::select(-chrm, -start, -end)

  feature_positions <- merged %>% dplyr::filter(!is.na(GENE))

  # split_to_codons <- function(x) {
  #   StartVec <- 0:(length(x)/3-1) * 3 + 1
  #   EndVec <- 1:(length(x)/3) * 3
  #   gene_seq <- x %>% unlist %>% paste0(collapse = "")
  #   codons <- stringr::str_sub(string = gene_seq, start = StartVec, end = EndVec)
  #   rep(codons, each = 3)
  # }

  feature_positions <- feature_positions %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("REF_CODON" = get_codons(REF_BASE))

  codons <- Biostrings::DNAStringSet(feature_positions$REF_CODON)
  aas <- Biostrings::translate(codons, no.init.codon = TRUE)
  aas <- as.character(aas, use.names = FALSE)

  feature_positions <- feature_positions %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("GENE_AA_POS" = rep(1:(dplyr::n()/3), each = 3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate("REF_AA" = aas)

  nonfeature_positions <- merged %>%
    dplyr::filter(is.na(GENE)) %>%
    dplyr::mutate("REF_CODON" = NA) %>%
    dplyr::mutate("REF_AA" = NA) %>%
    dplyr::mutate("GENE_AA_POS" = NA)

  all_ref <- dplyr::bind_rows(nonfeature_positions,
                       feature_positions) %>%
    dplyr::arrange(POS)

  all_ref$GENE[which(is.na(all_ref$GENE))] <- "non-genic"

  all_ref <- all_ref %>% tidyr::unite(col = "ref1",
                                      GENE, REF_AA,
                                      sep = "_",
                                      remove = FALSE) %>%
    tidyr::unite(col = "REF_IDENT",
          ref1, GENE_AA_POS,
          sep = "",
          remove = FALSE) %>%
    dplyr::select(-ref1) %>%
    dplyr::mutate("REF_IDENT" = stringr::str_replace_all(REF_IDENT, 
                                                         "non-genic_NANA", 
                                                         "non-genic"))

  all_ref <- all_ref %>%
    dplyr::select(POS,
                 REF_BASE,
                 GENE,
                 REF_CODON,
                 REF_AA,
                 GENE_AA_POS,
                 REF_IDENT) %>%
    dplyr::mutate("REF_BASE" = as.character(REF_BASE))

  all_ref
}
