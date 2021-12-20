#' Create MixVir reference genome object
#'
#' Uses a fasta genome and bed file defining features of interest (genes/ORFs) to create a data frame that's used as a reference to translate nucleotide data to amino acids and subsequently call variants/mutations from a sample.
#' @param genome fasta formatted genome file
#' @param feature.bed bed file defining features of interest (open reading frames to translate). Tab delimited with 4 columns (without column names):"chr", "start", "end", "feature_name".
#' @keywords reference
#' @export
#' @return A data frame with columns CHR,POS,REF_BASE,GENE,REF_CODON,REF_AA,GENE_AA_POS,REF_IDENT
#' @examples
#' create_ref()

create_ref <- function(genome, feature.bed) {

  features <- readr::read_tsv(feature.bed, col_names = FALSE)
  names(features) <- c("chrm", "start", "end", "GENE")

  #create df that includes row for every position in each feature/gene
  #allows for overlapping genes
  #currently assumes all genes are on the same strand
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

  #to analyze multiple chromosomes, need to merge on CHR+POS
  features <- features %>%
    tidyr::unite("chr_pos",
          chrm, POS,
          sep = "_")
  
  ref_genome <- ref_genome %>%
    tidyr::unite("chr_pos",
          CHR, POS,
          sep = "_",
          remove = FALSE)
  
  merged <- dplyr::left_join(x = ref_genome,
                             y = features,
                             by = "chr_pos") %>%
    dplyr::select(-start, -end, -chr_pos)

  #get positions associated with genes/ORFs and get codons they're associated with (each codon should be repeated 3 times)
  feature_positions <- merged %>% 
    dplyr::filter(!is.na(GENE)) %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("REF_CODON" = get_codons(REF_BASE))

  #translate codons to amino acids
  codons <- Biostrings::DNAStringSet(feature_positions$REF_CODON)
  aas <- Biostrings::translate(codons, no.init.codon = TRUE)
  aas <- as.character(aas, use.names = FALSE)

  #add columns with the relative amino acid position within the gene and the amino acid identity
  feature_positions <- feature_positions %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("GENE_AA_POS" = rep(1:(dplyr::n()/3), each = 3)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate("REF_AA" = aas)

  #create df with positions not associated with a gene/ORF
  nonfeature_positions <- merged %>%
    dplyr::filter(is.na(GENE)) %>%
    dplyr::mutate("REF_CODON" = NA) %>%
    dplyr::mutate("REF_AA" = NA) %>%
    dplyr::mutate("GENE_AA_POS" = NA)

  #merge the feature and non-feature positions together and sort by chr and position
  all_ref <- dplyr::bind_rows(nonfeature_positions,
                       feature_positions) %>%
    dplyr::arrange(CHR, POS)

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
    dplyr::select(CHR,
                  POS,
                 REF_BASE,
                 GENE,
                 REF_CODON,
                 REF_AA,
                 GENE_AA_POS,
                 REF_IDENT) %>%
    dplyr::mutate("REF_BASE" = as.character(REF_BASE))

  all_ref
}

