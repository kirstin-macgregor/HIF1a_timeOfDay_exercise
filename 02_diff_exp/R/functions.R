##########################################################################
## Project    : HIF1a KO clock ex 
## Data       : RNA-seq
## Author     : Kirstin MacGregor
## Date       : Project start: January 04, 2023. Last update: January 23, 2023
## Version    : v1.0
##########################################################################


## ------------------------------------------------------------------------
## General functions
## ------------------------------------------------------------------------

#' saving plots from list
#'
#' @param plot_list list of plots
#' @param save_path pathway for files to be saved to
#'
#' @return

plot_list_save_function <- function(plot_list, save_chr, save_path, width, height) {
  plot_names <- stringr::str_c(names(plot_list), save_chr, ".pdf")
  purrr::pwalk(list(plot_names, plot_list),
               ggplot2::ggsave,
               width = width,
               height = height,
               path = save_path
  )
}

#' save individual plot
#'
#' @param filename
#' @param plot
#' @param save_chr
#' @param save_path
#' @param width
#' @param height
#'
#' @return
plot_save_function <- function(filename, plot, save_path, width, height) {
  ggplot2::ggsave(
    filename,
    plot,
    width = width,
    height = height,
    path = save_path,
    bg = "white"
  )
}

##------------------------------------------------------------------------
## Data wrangling
##------------------------------------------------------------------------

#' Tidy DE output file for saving as xlsx
#'
#' @param data
#'
#' @return
tidy_DE_output_for_saving_fun <- function(data) {
  data |> dplyr::select(!c(comp, main_effect)) |> dplyr::arrange(FDR)
}

#' Tidy raw count data
#'
#' @param data
#'
#' @return tidy dataframe
tidy_raw_counts_fun <- function(data) {
  data |>
    tibble::column_to_rownames("Geneid") |>
    dplyr::select(!c("seqnames", "start", "end", "length", "strand"))
}


## ------------------------------------------------------------------------
## DESeq analysis
## ------------------------------------------------------------------------

#' annotate DESeq results object
#'
#' @param data
#'
#' @return annotated data frame of DESeq results
annotate_dat_fun <- function(data, comp_name, main_effect) {
  results <- as.data.frame(data) |>
    tibble::rownames_to_column("gene_id_version") |>
    dplyr::left_join(annotations, "gene_id_version") |>
    dplyr::rename(logFC = log2FoldChange, FDR = padj) |>
    dplyr::select(!c(gene_biotype, seq_coord_system, description, canonical_transcript)) |>
    dplyr::mutate(comp= paste(comp_name), main_effect= paste(main_effect))
}


## ------------------------------------------------------------------------
## DESeq2 model checks
## ------------------------------------------------------------------------

#' Histogram of p-values
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
pval_hist_fun <- function(data, title) {
  ggplot2::ggplot() +
    ggplot2::geom_histogram(data= data, ggplot2::aes(x= pvalue, fill= FDR > 0.1))+
    ggplot2::scale_fill_manual(values = c("#440154FF", "#808080", "#808080")) +
    ggplot2::facet_wrap(~comp)+
    ggplot2::theme_bw()+
    ggplot2::ggtitle(title)
}

