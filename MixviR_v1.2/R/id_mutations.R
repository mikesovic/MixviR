##removed AF in upstream script (consensus_summary.R) and replaced with ALT_COUNT
##made corresponding changes in this script.

#' Add Read Depths To Ref
#'
#' Add column for total read depths for a given sample to reference df
#' @param ref reference genome in "MixVir" format (genomic positions repeated for each associated feature they're associated with, etc.)
#' @param samp Character vector length 1 storing the name of a csv file with data for one sample. Should have columns CHROM,POS,REF,ALT,DP,AF. Others will be ignored.
#' @param samp.dir Directory storing sample file defined with parameter 'samp'.
#' @keywords depth
#' @return Data frame with cols "genomic_pos"	"ref_base"	"gene"	"ref_codon"	"ref_AA	GENE_AA_POS"	"ref_identity" "DP"
#' @export
#' @examples
#' add_depths_to_ref()

add_depths_to_ref <- function(ref, samp, samp.dir) {

  sample_data <- readr::read_csv(paste0(samp.dir, "/", samp))
  #has CHROM	POS	REF	ALT	DP	ALT_COUNT

  #get idx's where pos is duplicated - these will include indels
  dup_idx <- which(duplicated(sample_data$POS))
  dup_positions <- sample_data$POS[dup_idx]

  sample_data <- sample_data %>%
    dplyr::filter(!(POS %in% dup_positions & ALT == '.')) %>%
    dplyr::arrange(POS, desc(DP)) %>%
    dplyr::distinct(POS, .keep_all = TRUE) %>%
    dplyr::select(POS, DP) 
  
  ref <- dplyr::left_join(x = ref,
                      y = sample_data,
                      by = "POS")

  ref
}


#' ID SNV-based Amino Acid Changes
#'
#' Identify amino acid changes associated with single nucleotide variation. Changes associated with indels are identified in separate function. Used by call_mutations function.
#' @param variant.calls Data frame with cols POS, REF, ALT, AF (alt freq), DP (total read depth). Additional columns will be ignored.
#' @param ref reference genome in "MixVir" format (genomic positions repeated for each associated feature they're associated with, etc.)
#' @keywords snps
#' @return Data frame with cols "POS", "REF_BASE", "GENE", "REF_CODON", "REF_AA", "GENE_AA_POS", "REF_IDENT", "REF", "ALT", "ALT_freq", "ALT_COUNT", "samp_codon", "samp_AA", "samp_identity", "DP"
#' @export
#' @examples
#' id_snps()
id_snps <- function(variant.calls, ref) {
  variants <- variant.calls %>% dplyr::select(POS, REF, ALT, AF, ALT_COUNT)
  #names(variants) <- c("genomic_pos", "REF", "ALT", "ALT_freq", "ALT_COUNT")
  #cut indels down to first position - will deal with these separately in id_indels()
  variants$ALT <- stringr::str_sub(variants$ALT, start = 1L, end = 1L)

  #merge sample variants with reference on position and add in ref base where no variant
  all_samp <- dplyr::left_join(x = ref, y = variants, by = "POS")

  #add the reference bases to col 'ALT' where there is no variant
  na_idx <- which(is.na(all_samp$ALT))
  all_samp$ALT[na_idx] <- all_samp$REF_BASE[na_idx]

  #add column with codon each feature position is associated with
  #first get codons associated with annotated features
  features_w_codons <- all_samp %>% dplyr::filter(GENE != "non-genic") %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("samp_codon" = get_codons(ALT))

  #translate codons
  codons <- Biostrings::DNAStringSet(features_w_codons$samp_codon)
  aas <- Biostrings::translate(codons, no.init.codon = TRUE) %>% as.character(use.names = FALSE)
  features_w_codons$samp_AA <- aas
  features_w_codons <- features_w_codons %>%
    dplyr::mutate("samp_identity" = paste0(GENE,
                                           "_",
                                           samp_AA,
                                           GENE_AA_POS)) %>%
    dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                  REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                  samp_AA, samp_identity, DP)

  #get the positions not associated with features to bind back in
  nonfeature_positions <- all_samp %>%
    dplyr::filter(GENE == "non-genic") %>%
    dplyr::mutate("samp_codon" = rep(NA, dplyr::n())) %>%
    dplyr::mutate("samp_AA" = rep(NA, dplyr::n())) %>%
    dplyr::mutate("samp_identity" = paste0(GENE,
                                           "_",
                                           samp_AA,
                                           GENE_AA_POS)) %>%
    dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                  REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                  samp_AA, samp_identity, DP)

  rbind(nonfeature_positions, as.data.frame(features_w_codons))
}


#' ID Indel-based Amino Acid Changes
#'
#' Identify amino acid changes associaged with indel variation. Changes associated with SNVs are identified in separate function. Used by call_mutations function.
#' @param variant.calls Data frame with cols POS, REF, ALT, AF, DP. Additional columns will be ignored.
#' @param ref reference genome in "MixVir" format (genomic positions repeated for each associated feature they're associated with, etc.)
#' @keywords indel
#' @return Data frame with cols "POS", "REF_BASE", "GENE", "REF_CODON", "REF_AA", "GENE_AA_POS", "REF_IDENT", "REF", "ALT", "AF", "ALT_COUNT", "samp_codon", "samp_AA", "samp_identity", "DP"
#' @export
#' @examples
#' id_indels()

id_indels <- function(variant.calls, ref) {

  samp_calls_indels <- data.frame()

  variants <- variant.calls %>%
    dplyr::select(POS, REF, ALT, AF, ALT_COUNT)
  
  #get in-frame deletions and add them to 'samp_calls_indels' df
  dels_in_frame <- variants %>%
    dplyr::filter(stringr::str_length(REF) > 1) %>% #need to check on cases where both REF and ALT have lengths > 1
    dplyr::mutate("del_length" = stringr::str_length(REF)-1) %>%
    dplyr::mutate("aa_del_length" = del_length/3) %>%
    dplyr::filter(del_length %% 3 == 0)

  if (nrow(dels_in_frame) > 0) {
    dels_in_frame_adj <- dels_in_frame %>%
      dplyr::select(POS, REF, ALT, AF, ALT_COUNT, aa_del_length)

    dels_in_frame_w_ref <- dplyr::left_join(x = dels_in_frame_adj,
                                            y = ref,
                                            by = "POS")

    #Some in-frame deletions don't start at the first codon position. 
    #For the purpose of naming these deletions, the name will depend on where the "in-frame"
    #deletion starts relative to the first codon affected.
    #if the "in-frame" deletion starts at position 1 of the codon,
    #will use "GENE_AA_POS". If it starts at position 2 or 3, will use "GENE_AA_POS_adj"
    
    idx_to_adjust <- which(dels_in_frame_w_ref$codon_position != 3)
    dels_in_frame_w_ref$GENE_AA_POS[idx_to_adjust] <- dels_in_frame_w_ref$GENE_AA_POS[idx_to_adjust]+1
    
    GENE_AA_POSs <- dels_in_frame_w_ref$GENE_AA_POS
    del_gene_start_position <- dels_in_frame_w_ref$gene_base_num

    dels_in_frame_w_ref$ALT <- rep("del", nrow(dels_in_frame_w_ref))
    dels_in_frame_w_ref$samp_codon <- rep("del", nrow(dels_in_frame_w_ref))
    dels_in_frame_w_ref$samp_AA <- rep("del", nrow(dels_in_frame_w_ref))
    dels_in_frame_w_ref$samp_identity <- paste0("del",
                                                dels_in_frame_w_ref$GENE_AA_POS,
                                                "/",
                                                dels_in_frame_w_ref$GENE_AA_POS+dels_in_frame_w_ref$aa_del_length-1)

    del_starts <- dels_in_frame_w_ref$GENE_AA_POS
    del_ends <- dels_in_frame_w_ref$GENE_AA_POS+dels_in_frame_w_ref$aa_del_length-1
    del_name_edit_idx <- which(del_starts == del_ends)
    dels_in_frame_w_ref$samp_identity[del_name_edit_idx] <- gsub("/.+", "", dels_in_frame_w_ref$samp_identity[del_name_edit_idx])

    dels_in_frame_w_ref <- dels_in_frame_w_ref %>%
      dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)

    samp_calls_indels <- rbind(samp_calls_indels, dels_in_frame_w_ref)
  }

  #get frame-shift deletions and add them to 'samp_calls_indels' df
  dels_out_frame <- variants %>%
    dplyr::filter(stringr::str_length(REF) > 1) %>%
    dplyr::mutate("del_length" = stringr::str_length(REF)-1) %>%
    dplyr::mutate("aa_del_length" = del_length/3) %>%
    dplyr::filter(del_length %% 3 != 0)

  if (nrow(dels_out_frame) > 0) {
    dels_out_frame_adj <- dels_out_frame %>%
      dplyr::select(POS, REF, AF, ALT_COUNT, del_length)

    dels_out_frame_w_ref <- dplyr::left_join(x = dels_out_frame_adj,
                                             y = ref,
                                             by = "POS")

    dels_out_frame_w_ref$ALT <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_codon <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_AA <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_identity <- paste0("Fdel",
                                                 dels_out_frame_w_ref$GENE_AA_POS,
                                                 "/",
                                                 dels_out_frame_w_ref$del_length,
                                                 "bp")
    dels_out_frame_w_ref <- dels_out_frame_w_ref %>%
      dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)
    
    samp_calls_indels <- rbind(samp_calls_indels, dels_out_frame_w_ref)
  }

  #get in-frame insertions and add them to 'samp_calls_indels' df
  ins_in_frame <- variants %>%
    dplyr::filter(stringr::str_length(ALT) > 1) %>%
    dplyr::mutate("ins_length" = stringr::str_length(ALT)-1) %>%
    dplyr::mutate("aa_ins_length" = ins_length/3) %>%
    dplyr::filter(ins_length %% 3 == 0)

  if (nrow(ins_in_frame) > 0) {
    ins_in_frame_adj <- ins_in_frame %>%
      dplyr::select(POS, REF, AF, ALT_COUNT, aa_ins_length)

    ins_in_frame_w_ref <- dplyr::left_join(x = ins_in_frame_adj,
                                           y = ref,
                                           by = "POS")

    ins_in_frame_w_ref$ALT <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_codon <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_AA <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_identity <- paste0("ins",
                                               ins_in_frame_w_ref$GENE_AA_POS,
                                               "/",
                                               ins_in_frame_w_ref$GENE_AA_POS+ins_in_frame_w_ref$aa_ins_length-1)

    ins_in_frame_w_ref <- ins_in_frame_w_ref %>%
      dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)
    
    samp_calls_indels <- rbind(samp_calls_indels, ins_in_frame_w_ref)
  }

  #get frame-shift insertions and add them to 'samp_calls_indels' df
  ins_out_frame <- variants %>%
    dplyr::filter(stringr::str_length(ALT) > 1) %>%
    dplyr::mutate("ins_length" = stringr::str_length(ALT)-1) %>%
    dplyr::mutate("aa_ins_length" = ins_length/3) %>%
    dplyr::filter(ins_length %% 3 != 0)

  if (nrow(ins_out_frame) > 0) {
    ins_out_frame_adj <- ins_out_frame %>%
      dplyr::select(POS, REF, AF, ALT_COUNT, ins_length)

    ins_out_frame_w_ref <- dplyr::left_join(x = ins_out_frame_adj,
                                            y = ref,
                                            by = "POS")

    ins_out_frame_w_ref$ALT <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_codon <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_AA <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_identity <- paste0("Fins",
                                                ins_out_frame_w_ref$GENE_AA_POS,
                                                "/",
                                                ins_out_frame_w_ref$ins_length,
                                                "bp")
    ins_out_frame_w_ref <- ins_out_frame_w_ref %>%
      dplyr::select(POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)
    
    samp_calls_indels <- rbind(samp_calls_indels, ins_out_frame_w_ref)

  }

  samp_calls_indels
}

#' Identify Sample Amino Acid Changes
#'
#' Function to identify full set of amino acid changes in a given sample (includes changes based on both SNVs and indels)
#' @param sample.dir Path to directory with one tab-delimited file for each sample to analyze.
#' Each file should contain columns named POS, REF, ALT, AF, DP and is generally
#' a summary from a vcf file. Additional columns can be included and will be ignored.
#' POS = genomic position, REF = reference base, ALT = alternate base/allele, AF = alt frequency, DP = read depth at site
#' @param min.alt.freq Minimum frequency (0-1) for retaining alternate allele. Default = 0.01.
#' @param name.sep Character in sample names that separates the unique sample identifier (characters preceeding the separator) from any additional text. Only text preceeding the first instance of the character will be retained.
#' @param reference MixviR ref object. "Wuhan" uses pre-generated Sars-Cov2 ref genome. Otherwise, must provide arguments fasta.genome and bed to generate ref genome object on the fly.
#' @param write.mut.table Logical to indicated whether to write a text file that stores all mutations
#' for all samples analyzed to working directory. This is the same information contained in data frame
#' 'samp_mutations' created by this function. Default = FALSE
#' @param fasta.genome fasta formatted genome file (assumes only one chromosome)
#' @param bed bed file (4 columns, tab delimited, no column names: chr, start, end, feature_name) with info on features of interest (open reading frames to translate)
#' @keywords mutation
#' @return Data frame 'samp_mutations' containing amino acid changes observed for each sample.
#' @export
#' @examples
#' call_mutations()



#name.sep - character used to pull out unique identifier for each sample. The portion of the sample name preceeding the first instance of this character is retained, and should uniquely identify the sample.
call_mutations <- function(sample.dir,
                           min.alt.freq = 0.01, ###need to apply a filter before we call mutations. Otherwise, lots of noise gets included in calls and can affect real calls.
                           name.sep = "NULL",
                           reference = "Wuhan",
                           write.mut.table = FALSE,
                           fasta.genome,
                           bed) {

  samp_files <- dir(sample.dir)
  
  if ((reference == "Wuhan" | reference == "wuhan")) {
    ref_no_dp <- readr::read_tsv("https://raw.githubusercontent.com/mikesovic/IDI-AMSL/main/SC2_ref.tsv")
  }
  
  else {
    ref_no_dp <- create_ref(genome = fasta.genome, feature.bed = bed)
  }

  #df to store all mutation calls across all samples
  all_variants <- data.frame()

  for (file in samp_files) {
    curr_samp <- file

    if (!is.null(name.sep)) {
      curr_samp <- gsub(paste0("(.+?)", name.sep, "(.*)"), "\\1", file)
    }

    ref <- add_depths_to_ref(ref = ref_no_dp,
                             samp = file,
                             samp.dir = sample.dir)
    
    ref <- ref %>%
      dplyr::group_by(GENE) %>%
      dplyr::mutate("gene_base_num" = 1:dplyr::n()) %>%
      dplyr::mutate("codon_position" = dplyr::case_when(gene_base_num %% 3 == 0 ~ 3,
                                                        gene_base_num %% 3 == 1 ~ 1,
                                                        gene_base_num %% 3 == 2 ~ 2)) %>%
      dplyr::ungroup()
    
    all_variants_temp <- data.frame()
    print(curr_samp)

    sample_variants <- readr::read_csv(paste0(sample.dir, "/", file),
                                       col_types = cols()) %>%
      dplyr::filter(ALT != '.') %>%
      dplyr::select(POS, REF, ALT, ALT_COUNT, DP) %>%
      dplyr::mutate("AF" = ALT_COUNT/DP) %>%
      dplyr::filter(AF >= min.alt.freq)
      
    #determine if there are any positions with multiple mutations
    multiple_mutation_idx <- which(duplicated(sample_variants$POS))

    #if no sites with multiple mutations
    if (length(multiple_mutation_idx) == 0) {
      samp_calls_snv <- id_snps(variant.calls = sample_variants, ref = ref)

      #run function to identify indels and add them to 'all_variants' df
      samp_calls_indels <- id_indels(variant.calls = sample_variants, ref = ref)

      all_variants_temp <- rbind(all_variants_temp, samp_calls_snv, samp_calls_indels)
    } else {    #if one or more sites with multiple mutations

      dups_df <- sample_variants %>% dplyr::slice(multiple_mutation_idx)
      sample_variants <- sample_variants %>% dplyr::slice(-multiple_mutation_idx)
      #deal with sample_variants here the same way as above, but then
      #subsequently add in the dups_df.

      samp_calls_snv <- id_snps(variant.calls = sample_variants, ref = ref)

      #run function to identify indels and add them to 'all_variants' df
      samp_calls_indels <- id_indels(variant.calls = sample_variants, ref = ref)

      all_variants_temp <- rbind(all_variants_temp, samp_calls_snv, samp_calls_indels)
      
      #deal with multiple mutation sites
      while (nrow(dups_df) > 0) {
        dup_idx <- which(duplicated(dups_df$POS))
        if (length(dup_idx) > 0) {
          not_dups <- dups_df %>% dplyr::slice(-dup_idx) #work with these in current iteration
          dups_df <- dups_df %>% dplyr::slice(dup_idx)

          samp_calls_snv <- id_snps(variant.calls = not_dups, ref = ref)
          samp_calls_indels <- id_indels(variant.calls = not_dups, ref = ref)

          all_variants_temp <- rbind(all_variants_temp, samp_calls_snv, samp_calls_indels)

        } else {
          not_dups <- dups_df
          dups_df <- data.frame()

          samp_calls_snv <- id_snps(variant.calls = not_dups, ref = ref)
          samp_calls_indels <- id_indels(variant.calls = not_dups, ref = ref)

          all_variants_temp <- rbind(all_variants_temp, samp_calls_snv, samp_calls_indels)
        }
      }
    }

    #for the current sample, pull out mutations (differences in ref AA and sample AA)
    #add this to the master 'all_variants' df
    all_variants_temp <- all_variants_temp %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      dplyr::mutate("samp_identity" = stringr::str_replace_all(samp_identity, "non-genic_NANA", "non-genic")) %>%
      dplyr::filter(REF_IDENT != samp_identity) %>%
      dplyr::filter(!is.na(AF)) %>%
      dplyr::mutate("samp_name" = curr_samp) %>%
      dplyr::arrange(POS)

    all_variants <- rbind(all_variants, all_variants_temp)
  }

  #create a ALT_ID column to use for merging
  #have to do this separately for indels and nonindels b/c naming structures are different
  #then put them back together in 'all_variants'
  nonindels <- all_variants %>%
    dplyr::filter(ALT %in% c("A", "C", "T", "G", "stop", "Stop", "*", "STOP"))

  nonindels$ALT_ID <- paste0(nonindels$GENE,
                             "_",
                             nonindels$REF_AA,
                             nonindels$GENE_AA_POS,
                             nonindels$samp_AA)

  indels <- all_variants %>%
    dplyr::filter(ALT %in% c("del", "ins"))

  indels$ALT_ID <- paste0(indels$GENE,
                          "_",
                          indels$samp_identity)

  all_variants <- rbind(nonindels, indels)

  #clean up all_variants
  all_variants <- all_variants %>%
    dplyr::select(-REF_IDENT, -samp_identity) %>%
    dplyr::arrange(samp_name, POS)

  #create samp_mutations df and write it out - this is the final output
  samp_mutations <- all_variants %>%
    dplyr::select(samp_name, GENE, POS, ALT_ID, AF, ALT_COUNT, DP) %>%
    dplyr::rename("SAMP_NAME" = "samp_name")

  if (write.mut.table == TRUE) {
    write.table(samp_mutations, file = "sample_mutations.tsv",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }

  samp_mutations <<- samp_mutations
}
