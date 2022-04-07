library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

#get the reference df
reference <- read_tsv("https://raw.githubusercontent.com/mikesovic/IDI-AMSL/main/SC2_ref.tsv")

#get the called vcf for the sample
vcf <- read_tsv(paste0("vcfs/", args[1], ".vcf"),
                comment = "##") %>%
  select('#CHROM', POS, REF, ALT)

names(vcf) <- c("CHROM", "POS", "REF", "ALT")           

depths <- read_tsv(paste0("depths/", args[1], ".ldepth")) %>%
  select(SUM_DEPTH)

counts <- read_tsv(paste0("AD/", args[1], ".AD.FORMAT"), col_types = list("c", "i", "c"))
names(counts) <- c("CHROM", "POS", "ALLELE_COUNTS")
counts <- counts %>% select(ALLELE_COUNTS)

vcf_w_depth_and_cts <- cbind(vcf, depths, counts)

vcf_w_depth_and_freqs <- vcf_w_depth_and_cts %>%
  separate(col = ALLELE_COUNTS,
           into = c("REF_COUNT", "ALT_COUNT"),
           sep = ",",
           remove = FALSE) %>%
  select(CHROM, POS, REF, ALT, SUM_DEPTH, REF_COUNT, ALT_COUNT)

vcf_w_depth_and_freqs$REF_COUNT <- as.integer(vcf_w_depth_and_freqs$REF_COUNT)
vcf_w_depth_and_freqs$ALT_COUNT <- as.integer(vcf_w_depth_and_freqs$ALT_COUNT)

vcf_w_depth_and_freqs$ALT_COUNT <- replace_na(vcf_w_depth_and_freqs$ALT_COUNT, 0)

vcf_w_depth_and_freqs <- vcf_w_depth_and_freqs %>%
  select(CHROM, POS, REF, ALT, SUM_DEPTH, ALT_COUNT) %>%
  dplyr::rename("DP" = "SUM_DEPTH")

write.table(vcf_w_depth_and_freqs,
            file = paste0("consensus/", args[1], ".csv"),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
