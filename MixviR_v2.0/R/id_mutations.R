#' Convert VCF Samples to MixviR Input Format
#'
#' Create csv files with relevant contents of VCF to use as input for MixviR
#' @param vcf.dir path to directory containing one or more vcf files that will be converted to csv input format for MixviR call_mutations(). VCF's need to contain "DP" and "AD" flags in the FORMAT field. Should not contain any other files.
#' @param csv.path Path to directory where csv files will be written.
#' @keywords VCF
#' @export
#' @return csv file(s) with cols "CHR"	"POS"	"REF"	"ALT"	"DP"	"REF_COUNT" ALT_COUNT". Written to path defined with 'csv.path'.
#' @examples
#' bulk_vcf_to_mixvir()

bulk_vcf_to_mixvir <- function(vcf.dir = NULL, csv.path = NULL,
                               max.vcf.size = 1e+08){
  
    vcf_names <- dir(vcf.dir)
    
    for (infile in vcf_names) {
      
      vcf.dir <- gsub("/$", "", vcf.dir)
      
      #read in vcf
      vcf_obj <- vcfR::read.vcfR(paste0(vcf.dir, "/", infile), 
                                 limit = max.vcf.size,
                                 verbose = FALSE)
      
      #get df with total depths at each position
      depths <- vcfR::extract.gt(vcf_obj, "DP", as.numeric = TRUE) %>%
        as.data.frame()
      names(depths) <- c("DP")
      
      #get df with depths of each allele
      allele_depths <- vcfR::extract.gt(vcf_obj, "AD") %>%
        as.data.frame()
      names(allele_depths) <- c("ALLELE_COUNTS")
      allele_depths <- allele_depths %>%
        dplyr::mutate("ALLELE_COUNTS" = as.character(ALLELE_COUNTS)) %>%
        tidyr::separate(col = ALLELE_COUNTS,
                        into = c("REF_COUNT", "ALT_COUNT"),
                        sep = ",",
                        fill = "right") %>%
        dplyr::mutate("REF_COUNT" = as.integer(REF_COUNT)) %>%
        dplyr::mutate("ALT_COUNT" = as.integer(ALT_COUNT))
      allele_depths$ALT_COUNT <- tidyr::replace_na(allele_depths$ALT_COUNT, 0)
      
      
      #read in the data part of the vcf as a data frame
      consensus <- readr::read_tsv(paste0(vcf.dir, "/", infile),
                                   show_col_types = FALSE,
                                   comment = "##") %>%
        dplyr::select('#CHROM', POS, REF, ALT) %>%
        dplyr::rename("CHR" = '#CHROM')
      
      #add in the total depths/allele depths - this gives the original "consensus" format for MixviR
      consensus <- cbind(consensus, depths, allele_depths)
      
      csv.path <- gsub("/$", "", csv.path)
      infile <- gsub(".vcf$", "", infile)
      write.table(consensus, file = paste0(csv.path, "/", infile, ".csv"),
                  sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    }
  }

  
#' Convert Sample VCF to MixviR Input Format
#'
#' Create data frame with relevant contents of VCF
#' @param infile path to a vcf file
#' @param max.vcf.size Max memory usage (in bytes) allowed when reading in vcf file (from vcfR) 
#' @keywords VCF
#' @return Data frame with cols "CHR"	"POS"	"REF"	"ALT"	"DP" "REF_COUNT"	"ALT_COUNT"
#' @examples
#' vcf_to_mixvir()

vcf_to_mixvir <- function(infile, max.vcf.size = 1e+08) {
  
    #read in vcf
    vcf_obj <- vcfR::read.vcfR(infile, 
                           limit = max.vcf.size,
                           verbose = FALSE)
    
    #get df with total depths at each position
    depths <- vcfR::extract.gt(vcf_obj, "DP", as.numeric = TRUE) %>%
      as.data.frame()
    names(depths) <- c("DP")
    
    #get df with depths of each allele
    allele_depths <- vcfR::extract.gt(vcf_obj, "AD") %>%
      as.data.frame()
    names(allele_depths) <- c("ALLELE_COUNTS")
    allele_depths <- allele_depths %>%
      dplyr::mutate("ALLELE_COUNTS" = as.character(ALLELE_COUNTS)) %>%
      tidyr::separate(col = ALLELE_COUNTS,
                      into = c("REF_COUNT", "ALT_COUNT"),
                      sep = ",",
                      fill = "right") %>%
      dplyr::mutate("REF_COUNT" = as.integer(REF_COUNT)) %>%
      dplyr::mutate("ALT_COUNT" = as.integer(ALT_COUNT))
    allele_depths$ALT_COUNT <- tidyr::replace_na(allele_depths$ALT_COUNT, 0)
    
    
    #read in the data part of the vcf as a data frame
    consensus <- readr::read_tsv(infile,
                    comment = "##",
                    show_col_types = FALSE) %>%
      dplyr::select('#CHROM', POS, REF, ALT) %>%
      dplyr::rename("CHR" = '#CHROM')
    
    #add in the total depths/allele depths - this gives the original "consensus" format for MixviR
    consensus <- cbind(consensus, depths, allele_depths)
    
    return(consensus)
}


#' Add Read Depths To Ref
#'
#' Add column of total read depths for a given sample to reference df
#' @param ref reference genome in "MixVir" format (from create_ref() function)
#' @param samp.variants Data frame produced by function vcf_to_mixvir(). Contains columns "CHR", "POS", "REF", "ALT", "DP", "ALT_COUNT".
#' @keywords depth
#' @return Original 'ref' data frame with depth (DP) column added: cols "genomic_pos"	"ref_base"	"gene"	"ref_codon"	"ref_AA"	"GENE_AA_POS"	"ref_identity" "DP"
#' @examples
#' add_depths_to_ref()

add_depths_to_ref <- function(ref, samp.variants) {

  #use combination of chr+pos to deal with multiple chromosomes
  sample_data <- samp.variants %>%
    dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))
    
  #get idx's where pos is duplicated - these will include indels
  dup_idx <- which(duplicated(sample_data$chr_pos))
  dup_positions <- sample_data$chr_pos[dup_idx]

  sample_data <- sample_data %>%
    dplyr::filter(!(chr_pos %in% dup_positions & ALT == '.')) %>%
    dplyr::arrange(CHR, POS, desc(DP)) %>%
    dplyr::distinct(chr_pos, .keep_all = TRUE) %>%
    dplyr::select(chr_pos, DP) 
  
  ref <- ref %>%
    dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))
  
  ref <- dplyr::left_join(x = ref,
                      y = sample_data,
                      by = "chr_pos") %>%
    dplyr::select(-chr_pos)

  return(ref)
}


#' ID SNV-based Amino Acid Changes
#'
#' Identify amino acid changes associated with single nucleotide variation. Changes associated with indels are identified in separate function. Used by call_mutations() function.
#' @param variant.calls Data frame with cols POS, REF, ALT, AF (alt freq), DP (total read depth).
#' @param ref reference genome in "MixVir" format (from create_ref() function)
#' @keywords snps
#' @return Data frame that includes amino acid calls based on SNP/SNV variants observed in sample. Contains cols "POS", "REF_BASE", "GENE", "REF_CODON", "REF_AA", "GENE_AA_POS", "REF_IDENT", "REF", "ALT", "ALT_freq", "ALT_COUNT", "samp_codon", "samp_AA", "samp_identity", "DP"
#' @examples
#' id_snps()
id_snps <- function(variant.calls, ref) {
  variants <- variant.calls %>% 
    dplyr::select(CHR, POS, REF, ALT, AF, ALT_COUNT) %>%
    dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))
  
  #cut indels down to first position - will deal with these separately in id_indels()
  variants$ALT <- stringr::str_sub(variants$ALT, start = 1L, end = 1L)

  ref <- ref %>%
    dplyr::mutate("chr_pos" = paste0(CHR, "_", POS)) %>%
    dplyr::select(-CHR, -POS)
  #merge sample variants with reference on chr+position and add in ref base where no variant
  all_samp <- dplyr::left_join(x = ref, y = variants, by = "chr_pos")

  #add the reference bases to col 'ALT' where there is no variant
  na_idx <- which(is.na(all_samp$ALT))
  all_samp$ALT[na_idx] <- all_samp$REF_BASE[na_idx]

  #add column with codon each feature position is associated with
  #first get codons associated with annotated features
  features_w_codons <- all_samp %>% 
    dplyr::filter(GENE != "non-genic") %>%
    dplyr::group_by(GENE) %>%
    dplyr::mutate("samp_codon" = get_codons(ALT))

  #translate codons
  codons <- Biostrings::DNAStringSet(features_w_codons$samp_codon)
  aas <- Biostrings::translate(codons, no.init.codon = TRUE) %>% 
    as.character(use.names = FALSE)
  features_w_codons$samp_AA <- aas
  features_w_codons <- features_w_codons %>%
    dplyr::mutate("samp_identity" = paste0(GENE,
                                           "_",
                                           samp_AA,
                                           GENE_AA_POS)) %>%
    dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
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
    dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
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
#' @return Data frame that includes amino acid calls based on indel variants observed in sample. Contains cols "POS", "REF_BASE", "GENE", "REF_CODON", "REF_AA", "GENE_AA_POS", "REF_IDENT", "REF", "ALT", "ALT_freq", "ALT_COUNT", "samp_codon", "samp_AA", "samp_identity", "DP"
#' @examples
#' id_indels()

id_indels <- function(variant.calls, ref) {

  samp_calls_indels <- data.frame()

  variants <- variant.calls %>%
    dplyr::select(CHR, POS, REF, ALT, AF, ALT_COUNT)
  
  #get in-frame deletions and add them to 'samp_calls_indels' df
  dels_in_frame <- variants %>%
    dplyr::filter(stringr::str_length(REF) > 1) %>% #need to check on cases where both REF and ALT have lengths > 1
    dplyr::mutate("del_length" = stringr::str_length(REF)-1) %>%
    dplyr::mutate("aa_del_length" = del_length/3) %>%
    dplyr::filter(del_length %% 3 == 0)  #select for just in-frame deletions

  ref <- ref %>%
    dplyr::mutate("chr_pos" = paste0(CHR, "_", POS)) %>%
    dplyr::select(-CHR, -POS)
  
  if (nrow(dels_in_frame) > 0) {
    dels_in_frame_adj <- dels_in_frame %>%
      dplyr::select(CHR, POS, REF, ALT, AF, ALT_COUNT, aa_del_length) %>%
      dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))
    
    
    dels_in_frame_w_ref <- dplyr::left_join(x = dels_in_frame_adj,
                                            y = ref,
                                            by = "chr_pos")

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

    #clean up names of 1-amino acid deletions - i.e. from del50/50 to del50
    del_starts <- dels_in_frame_w_ref$GENE_AA_POS
    del_ends <- dels_in_frame_w_ref$GENE_AA_POS+dels_in_frame_w_ref$aa_del_length-1
    del_name_edit_idx <- which(del_starts == del_ends)
    dels_in_frame_w_ref$samp_identity[del_name_edit_idx] <- gsub("/.+", "", dels_in_frame_w_ref$samp_identity[del_name_edit_idx])

    dels_in_frame_w_ref <- dels_in_frame_w_ref %>%
      dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)

    samp_calls_indels <- rbind(samp_calls_indels, dels_in_frame_w_ref)
  }

  #get frame-shift deletions and add them to 'samp_calls_indels' df
  dels_out_frame <- variants %>%
    dplyr::filter(stringr::str_length(REF) > 1) %>%
    dplyr::mutate("del_length" = stringr::str_length(REF)-1) %>%
    dplyr::mutate("aa_del_length" = del_length/3) %>%
    dplyr::filter(del_length %% 3 != 0) #select any out-of-frame deletions

  if (nrow(dels_out_frame) > 0) {
    dels_out_frame_adj <- dels_out_frame %>%
      dplyr::select(CHR, POS, REF, AF, ALT_COUNT, del_length) %>%
      dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))

    dels_out_frame_w_ref <- dplyr::left_join(x = dels_out_frame_adj,
                                             y = ref,
                                             by = "chr_pos")

    dels_out_frame_w_ref$ALT <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_codon <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_AA <- rep("del", nrow(dels_out_frame_w_ref))
    dels_out_frame_w_ref$samp_identity <- paste0("Fdel",
                                                 dels_out_frame_w_ref$GENE_AA_POS,
                                                 "/",
                                                 dels_out_frame_w_ref$del_length,
                                                 "bp")
    
    dels_out_frame_w_ref <- dels_out_frame_w_ref %>%
      dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
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
      dplyr::select(CHR, POS, REF, AF, ALT_COUNT, aa_ins_length) %>%
      dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))

    ins_in_frame_w_ref <- dplyr::left_join(x = ins_in_frame_adj,
                                           y = ref,
                                           by = "chr_pos")

    ins_in_frame_w_ref$ALT <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_codon <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_AA <- rep("ins", nrow(ins_in_frame_w_ref))
    ins_in_frame_w_ref$samp_identity <- paste0("ins",
                                               ins_in_frame_w_ref$GENE_AA_POS,
                                               "/",
                                               ins_in_frame_w_ref$GENE_AA_POS+ins_in_frame_w_ref$aa_ins_length-1)

    ins_in_frame_w_ref <- ins_in_frame_w_ref %>%
      dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
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
      dplyr::select(CHR, POS, REF, AF, ALT_COUNT, ins_length) %>%
      dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))

    ins_out_frame_w_ref <- dplyr::left_join(x = ins_out_frame_adj,
                                            y = ref,
                                            by = "chr_pos")

    ins_out_frame_w_ref$ALT <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_codon <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_AA <- rep("ins", nrow(ins_out_frame_w_ref))
    ins_out_frame_w_ref$samp_identity <- paste0("Fins",
                                                ins_out_frame_w_ref$GENE_AA_POS,
                                                "/",
                                                ins_out_frame_w_ref$ins_length,
                                                "bp")
    ins_out_frame_w_ref <- ins_out_frame_w_ref %>%
      dplyr::select(CHR, POS, REF_BASE, GENE, REF_CODON, REF_AA, GENE_AA_POS,
                    REF_IDENT, REF, ALT, AF, ALT_COUNT, samp_codon,
                    samp_AA, samp_identity, DP)
    
    samp_calls_indels <- rbind(samp_calls_indels, ins_out_frame_w_ref)

  }

  samp_calls_indels
}

#' Identify Amino Acid Changes From A Potentially Mixed Sample
#'
#' Identify full set of amino acid changes from one or more samples (includes changes based on both SNVs and indels)
#' @param sample.dir Path to directory containing vcf or csv files for each sample to analyze. VCF's need to contain "DP" and "AD" flags in the FORMAT field. Should not contain any other files.
#' @param csv.input Logical to indicate whether files in *sample.dir* directory are in vcf or csv format. All files must be of the same format. If csv, they must contain columns named: "CHR"	"POS"	"REF"	"ALT"	"DP"	"ALT_COUNT". See the `batch_vcf_to_mixvir()`` function to convert vcfs to csv format (csv may be slightly faster for running call_mutations()). Default is FALSE (input is in vcf format).
#' @param min.alt.freq Minimum frequency (0-1) for retaining alternate allele. Default = 0.01. Extremely low values (i.e. zero) are not recommended here - see vignette for details.
#' @param name.sep Character in input file names that separates the unique sample identifier (characters preceeding the separator) from any additional text. Only text preceeding the first instance of the character will be retained and used as the sample name.
#' @param reference MixviR ref object (created with `create_ref()``). "Wuhan" uses pre-generated Sars-Cov2 ref genome. Otherwise, must provide arguments *fasta.genome* and *bed* to generate ref genome object on the fly.
#' @param write.mut.table Logical indicating whether to write the 'samp_mutations' data frame (see "Value" below) to a text file in the current working directory. Default = FALSE.
#' @param fasta.genome fasta formatted reference genome file.
#' @param bed bed file defining features of interest (open reading frames to translate). Should be tab delimited and have 4 columns (no column names): chr, start, end, feature_name
#' @keywords mutation
#' @return Object 'samp_mutations' stored in the global environment. Specifically, a data frame containing amino acid changes observed for each sample, positions of the underlying mutations, and other information.
#' @export
#' @examples
#' call_mutations()

call_mutations <- function(sample.dir = NULL,
                           csv.infiles = FALSE,
                           min.alt.freq = 0.01, ###need to apply a filter before we call mutations. Otherwise, lots of noise gets included and can affect real calls.
                           name.sep = NULL,
                           reference = "custom",
                           write.mut.table = FALSE,
                           fasta.genome,
                           bed) {

  
  samp_files <- dir(sample.dir)
  
  if ( (reference !="Wuhan") && ((is.null(fasta.genome) || is.null(bed))) ) {
    print("Error: call_mutations() needs a reference genome and associated annotations. Must provide fasta.genome and bed files, or set reference option to Wuhan")
  }
  
  if (reference == "Wuhan") {
    ref_no_dp <- readr::read_tsv("https://raw.githubusercontent.com/mikesovic/MixviR/main/reference_files/SC2_wuhan.tsv", show_col_types = FALSE)
    print("Using Wuhan reference")
  } else {
     ref_no_dp <- create_ref(genome = fasta.genome, feature.bed = bed)
  }

  #df to store all mutation calls across all samples
  all_variants <- data.frame()

  for (file in samp_files) {
    curr_samp <- file

    if (!is.null(name.sep)) {
      curr_samp <- gsub(paste0("(.+?)", name.sep, "(.*)"), "\\1", file)
    }

    #trim off trailing forward slash if it's in the sample.dir path
    sample.dir <- gsub("/$", "", sample.dir)
    
    if (csv.infiles == FALSE) {
      #file is currently a vcf - need to get it in to "MixviR/csv format"
      variants_df <- vcf_to_mixvir(infile = paste0(sample.dir, "/", file))
    } else{
      variants_df <- readr::read_csv(file = paste0(sample.dir, "/", file, show_col_types = FALSE))
    }
    
    #add coverages at each position for the current sample to the ref object/data frame
    ref <- add_depths_to_ref(ref = ref_no_dp,
                             samp.variants = variants_df)
    
    #get the relative codon position for each base in each gene
    #i.e. does a given base sit in the 1st, 2nd, or 3rd position of its codon?
    #this is used later for naming indels
    ref <- ref %>%
      dplyr::group_by(GENE) %>%
      dplyr::mutate("gene_base_num" = 1:dplyr::n()) %>%
      dplyr::mutate("codon_position" = dplyr::case_when(gene_base_num %% 3 == 0 ~ 3,
                                                        gene_base_num %% 3 == 1 ~ 1,
                                                        gene_base_num %% 3 == 2 ~ 2)) %>%
      dplyr::ungroup()
    
    all_variants_temp <- data.frame()
    print(paste0("Calling mutations: ", curr_samp))
    
    sample_variants <- variants_df %>%
      dplyr::filter(ALT != '.') %>%
      dplyr::select(CHR, POS, REF, ALT, ALT_COUNT, DP) %>%
      dplyr::mutate("AF" = ALT_COUNT/DP) %>%
      dplyr::filter(AF >= min.alt.freq) %>%
      dplyr::mutate("chr_pos" = paste0(CHR, "_", POS))
      
    #determine if there are any positions with multiple mutations
    multiple_mutation_idx <- which(duplicated(sample_variants$chr_pos))

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
        dup_idx <- which(duplicated(dups_df$chr_pos))
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
      dplyr::arrange(CHR, POS)

    all_variants <- rbind(all_variants, all_variants_temp)
  }

  #create an ALT_ID column to use for merging
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
    dplyr::arrange(samp_name, CHR, POS)

  #create samp_mutations df and write it out - this is the final output
  samp_mutations <- all_variants %>%
    dplyr::select(samp_name, CHR, POS, GENE, ALT_ID, AF, ALT_COUNT, DP) %>%
    dplyr::rename("SAMP_NAME" = "samp_name")

  if (write.mut.table == TRUE) {
    write.table(samp_mutations, file = "sample_mutations.tsv",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
  }

  samp_mutations <<- samp_mutations
}

