#get list of all sample names
#get list of all odh mutations of interest
#edit deletions in samp_mutations to make them individual where they aren't
#copy samp_mutations df and clean - call samp_muts_odh
  #need cols 'samp_name', 'ALT_ID', 'ALT_COUNT', 'TOTAL_depth'
  #create mergable col 'merge_col' = samp_name+ALT_ID
  #delete samp_name and ALT_ID

#loop over samples
#for each sample...
  #filter samp_muts_odh for current sample -> samp_muts_odh_filtered
  #create 'merge_col' current_samp_name+all_mutation_IDs in new 'all_muts_df'
  #left join all_muts_df and filtered samp_muts_odh_filtered

#' Write Table In New ODH Format
#'
#' Write a table with read counts for each position
#' @param sample.dir Path to directory with one tab-delimited file for each sample to analyze.
#' Each file should contain columns named POS, REF, ALT, AF, DP and is generally
#' a summary from a vcf file. Additional columns can be included and will be ignored.
#' POS = genomic position, REF = reference base, ALT = alternate base/allele, AF = alt frequency, DP = read depth at site
#' @param mut.positions Data frame that defines amino acid changes of interest and associated mutation positions
#' @param name.sep Character in sample names that separates the unique sample identifier (characters preceeding the separator) from any additional text. Only text preceeding the first instance of the character will be retained.
#' @param reference Reference genome information in MixVir format.
#' @param site.meta Comma separated table with metadata for each sample site. Cols codes, Site_ID, Facility_Name, LAB, replicate
#' Should have two columns - 'MUT' and 'POS'
#' @keywords ODH
#' @export
#' @examples
#' write_ODH_table()

write_ODH_table <- function(sample.dir,
                            mut.positions = "https://raw.githubusercontent.com/mikesovic/IDI-AMSL/main/ODH_mut_pos_master_split_dels_6-3-21.txt",
                            name.sep = "NULL",
                            reference = "https://raw.githubusercontent.com/mikesovic/IDI-AMSL/main/SC2_ref.tsv",
                            site.meta = "https://raw.githubusercontent.com/mikesovic/IDI-AMSL/main/OSU_SequencingSites_CODEX.csv") {

  samp_files <- dir(sample.dir)
  print(samp_files)
  #get everything except deletions of length > 1
  samp_muts_odh <- samp_mutations %>%
    dplyr::select(samp_name, ALT_ID, ALT_COUNT, TOTAL_depth) %>%
    dplyr::slice(-grep("_del.+/.+", samp_mutations$ALT_ID))

  if (length(grep("_del.+/.+", samp_mutations$ALT_ID)) > 0) {
    #get deletions of length >1 and split them
    deletions <- samp_mutations %>%
      dplyr::slice(grep("_del.+/.+", samp_mutations$ALT_ID)) %>%
      dplyr::filter(gene != "non-genic") %>%
      dplyr::select(samp_name, ALT_ID, ALT_COUNT, TOTAL_depth) %>%
      tidyr::separate(col = ALT_ID,
                      into = c("gene", "ID"),
                      sep = "_") %>%
      dplyr::mutate("ID" = gsub("del", "", ID)) %>%
      tidyr::separate(col = ID,
                      into = c("start_aa", "end_aa"),
                      sep = "/")

    #think deletions with sizes that are multiples of 3 that are in non-genic regions
    #are getting lost in the above, so adding them back in.
    nongenic_check <- samp_mutations %>%
      dplyr::slice(grep("_del.+/.+", samp_mutations$ALT_ID)) %>%
      dplyr::filter(gene == "non-genic") %>%
      dplyr::select(samp_name, ALT_ID, ALT_COUNT, TOTAL_depth)

    if (nrow(nongenic_check) > 0) {
      samp_muts_odh <- dplyr::bind_rows(samp_muts_odh, nongenic_check)
    }

    aa_positions <- mapply(FUN = function(start, end) {
      start:end },
      start = deletions$start_aa,
      end = deletions$end_aa,
      SIMPLIFY = FALSE
    )
    names(aa_positions) <- NULL

    collapse_fun <- function(x) {
      x <- unlist(x)
      paste(x, collapse = ",")
    }

    aa_positions_collapsed <- lapply(aa_positions, collapse_fun)
    deletions$del <- aa_positions_collapsed

    deletions <- deletions %>%
      tidyr::separate_rows(del, sep = ",") %>%
      dplyr::mutate("ALT_ID" = paste0(gene, "_del", del)) %>%
      dplyr::select(samp_name, ALT_ID, ALT_COUNT, TOTAL_depth)

    samp_muts_split_indel <- dplyr::bind_rows(samp_muts_odh, deletions)
  } else {
    samp_muts_split_indel <- samp_muts_odh
  }

  odh_mutations <- readr::read_tsv(mut.positions) %>%
    dplyr::rename("ALT_ID" = "MUT",
                  "genomic_pos" = "POS")

  #this stores all mutation calls across all samples
  all_variants <- data.frame()


  for (file in samp_files) {

    curr_samp <- file

    if (!is.null(name.sep)) {
      curr_samp <- gsub(paste0("(.+?)", name.sep, "(.*)"), "\\1", file)
    }

    samp_muts_odh_filt <- samp_muts_split_indel %>%
      dplyr::filter(samp_name == curr_samp)

    #get the set of mutations of interest that didn't appear in the sample
    not_in_samp <- dplyr::anti_join(x = odh_mutations,
                                    y = samp_muts_odh_filt,
                                    by  = "ALT_ID")

    ref_df <- readr::read_tsv(reference)

    #get depths for positions associated with each mutation of interest not observed in the current sample
    ref_w_depth <- add_depths_to_ref(ref = ref_df,
                                     samp = file,
                                     samp.dir = sample.dir)

    depths <- ref_w_depth %>%
      dplyr::select(genomic_pos, DP)

    not_in_samp <- dplyr::left_join(x = not_in_samp,
                                    y = depths,
                                    by = "genomic_pos")

    not_in_samp <- not_in_samp %>%
      dplyr::mutate("samp_name" = curr_samp,
                    "ALT_COUNT" = 0,
                    "TOTAL_depth" = DP) %>%
      dplyr::select(samp_name, ALT_ID, ALT_COUNT, TOTAL_depth)

    all_variants <- dplyr::bind_rows(all_variants, samp_muts_odh_filt, not_in_samp)
  }

  all_variants_summary <- all_variants %>%
    dplyr::group_by(samp_name, ALT_ID) %>%
    dplyr::summarise("alt_reads" = mean(ALT_COUNT),
                     "total_reads" = mean(TOTAL_depth)) %>%
    tidyr::separate(col = "ALT_ID",
                    into = c("gene", "mutation"),
                    sep = "_")

  all_variants_summary$total_reads <- tidyr::replace_na(all_variants_summary$total_reads, 0)

  all_variants_summary$mutation <- gsub("(F?del)(.+)", "\\2\\1", all_variants_summary$mutation)
  all_variants_summary$mutation <- gsub("(F?ins)(.+)", "\\2\\1", all_variants_summary$mutation)

  #all_variants_summary$laboratory <- "OSU_AMSL"
  site_df <- readr::read_csv(site.meta, col_names = TRUE)

  all_variants_summary <- all_variants_summary %>%
    tidyr::separate(col = samp_name,
                    into = c("codex", "date"),
                    sep = "-")

  all_variants_summary <- dplyr::left_join(x = all_variants_summary,
                                           y = site_df,
                                           by = "codex") %>%
    dplyr::select(LAB, Site_ID, date, Facility_Name, replicate, gene, mutation, alt_reads, total_reads) %>%
    tidyr::unite(col = "sample_id",
                 Site_ID, date,
                 sep = "_",
                 remove = FALSE) %>%
    dplyr::select(LAB, sample_id, Facility_Name, replicate, gene, mutation, alt_reads, total_reads)

  write.table(all_variants_summary,
              file = "odh_mutation_read_counts.csv",
              sep = ",",
              row.names = FALSE,
              quote = FALSE)

}
