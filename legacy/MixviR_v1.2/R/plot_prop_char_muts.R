#' Plot Freqencies Of Lineage-Characteristic Mutations
#'
#' For each variant/lineage+sample combination, generates a barplot representing frequencies of each amino acid change associated with the lineage/variant.
#' @param genes Character vector defining genes to include in analysis. Default = "all".
#' @param min.read.depth Minimum number of reads required at the site to include (currently total reads)
#' @param lineage.def Tab-delimited file containing an a priori set of amino acid changes associated with variants/lineages of interest. Amino acid changes in the sample are matched to these. Choose from 'CDC' (default), 'Outbreak', 'Outbreak_May2021', or define your own.
#' @param variant.thresh Value (0-1) that defines horizontal reference line on plots
#' @param grid.dim Defines number of rows and columns plots for individual samples will be divided into for a given variant/lineage. Default: "c(1,1)" gives one plot per page.
#' @param lineage.def.name Name to use in plot titles to define the source of the lineage/variant definitions. Defaults to name for lineage.def argument.
#' @keywords plot
#' @import ggplot2
#' @export
#' @examples
#' plot_prop_charac_muts()

plot_prop_charac_muts <- function(genes = "all",
                                  min.read.depth = 50,
                                  lineage.def = "CDC",
                                  variant.thresh = 0.7,
                                  grid.dim = c(1,1),
                                  lineage.def.name = NULL) {

  if (is.null(lineage.def.name)) {
    lineage.def.name <- lineage.def
  }

  if (sum(grid.dim) > 2) {
    product_check <- grid.dim[1]*grid.dim[2]
  }

  #they want to plot by lineage (1 lineage per plot), so we'll loop over this -
  #it should be all the possible lineages (not just lineages present)
  all_uniq_lineages <- get.uniq.lineages(lineage.def = lineage.def)

  characterized_muts_samp <- merge_samp_muts_w_characteristics(lineage.def = lineage.def)

  all_uniq_samps <- characterized_muts_samp %>%
    dplyr::distinct(sample) %>%
    unlist()

  #if a gridded plot, check to make sure there are enough spaces for the number of samples
  if (sum(grid.dim) > 2) {
    product_check <- grid.dim[1]*grid.dim[2]
    if (product_check < length(all_uniq_samps)) {
      stop("Trying to create gridded plots, but (num rows * num cols) is less than number of samples to plot.")
    }
  }

  #summarize by lineage+samp+date
  #create the summary with counts of mutations w/RD > min.read.depth, total possible mutations, and proportion

  characterized_muts_samp$read_depth <- tidyr::replace_na(characterized_muts_samp$read_depth, 0)

  summary <- characterized_muts_samp %>%
    dplyr::mutate("read_depth_high" = ifelse(read_depth > min.read.depth, 1, 0)) %>%
    dplyr::group_by(Lineage, sample, date) %>%
    dplyr::summarise("n_characterized_muts" = dplyr::n(),
              "muts_in_samp" = sum(read_depth_high),
              "proportion" = muts_in_samp/n_characterized_muts)

  pdf(file = paste0(lineage.def, "_prop_charactistic_mutations.pdf"))

  for (lineage in all_uniq_lineages) {
    samp_counter <- 0
    plots_list <- vector('list', length(all_uniq_samps))

    for (samp in all_uniq_samps) {
      samp_counter <- samp_counter + 1

      data <- summary %>%
        dplyr::filter(Lineage == lineage & sample == samp)


      if (sum(grid.dim) > 2) {
        plots_list[[samp_counter]] <- local({

          samp_counter <- samp_counter

          plot <- data %>% ggplot2::ggplot(aes(x = date, y = proportion)) +
            geom_col() +
            ylim(c(0,1)) +
            xlab(NULL) +
            ylab("Proportion") +
            ggtitle(samp) +
            #subtitle = "Proportion Of Lineage-Characteristic Mutations W/Read Depth > 50") +
            theme(axis.text.x = element_text(angle = 45,
                                             vjust = 0.9,
                                             hjust = 1,
                                             face = "bold"),
                  axis.text.y = element_text(size = 6),
                  axis.title.y = element_text(size = 6),
                  plot.title = element_text(size=10)) +
            geom_hline(aes(yintercept = variant.thresh), color = "grey")

        })
      }  else {
        plot <- data %>% ggplot2::ggplot(aes(x = date, y = proportion)) +
          geom_col() +
          ylim(c(0,1)) +
          xlab(NULL) +
          ylab("Proportion") +
          ggtitle(paste0(lineage, " : ", samp),
                  subtitle = "Proportion Of Lineage-Characteristic Mutations W/Read Depth > 50") +
          theme(axis.text.x = element_text(angle = 45,
                                           vjust = 0.9,
                                           hjust = 1,
                                           face = "bold")) +
          geom_hline(aes(yintercept = variant.thresh), color = "grey")

        print(plot)
      }
    }

    if (sum(grid.dim) > 2) {
      p <- patchwork::wrap_plots(plots_list, nrow = grid.dim[1], ncol = grid.dim[2]) +
        patchwork::plot_annotation(title = paste0("Proportion ", lineage, "-Like Mutations (", lineage.def.name, ")"),
                        subtitle = paste0("Reference line at ", variant.thresh*100, "% indicates threshold for presence of lineage"))

      print(p)
    }
  }
  dev.off()
}
