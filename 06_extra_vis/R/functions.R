##########################################################################
## Project    : Functions for HIF1a analysis
## Data       : transcriptomics
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

##------------------------------------------------------------------------
## GO term stuff
##------------------------------------------------------------------------

#' Return list of genes within a specified GO term.
#'
#' @param GO_term
#'
#' @return
GO_ToI_fun <- function(GO_term) {
  results <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = c(GO_term), columns = c("ENSEMBL"), keytype = "GOALL")
  gene_symbols <- unique(results$ENSEMBL)
}

#' Count number of DE genes in a specified GO term per group
#'
#' @param data
#' @param GO_term
#'
#' @return
count_DE_genes_in_ToI <- function(data, GO_term) {
  as.data.frame(data) |>
    dplyr::filter(
      gene_id %in% GO_ToI_fun(GO_term),
      FDR < 0.1
    ) |>
    dplyr::group_by(comp, .drop = FALSE) |>
    dplyr::summarise(n = length(gene_id))
}

#' Bar plot for number of DE genes in a specified GO term per group
#'
#' @param data
#' @param title
#' @param comp
plot_DE_genes_in_ToI <- function(data, title, comp, ylim) {
  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_col(ggplot2::aes(x = comp, y = n, fill = comp)) +
    ggplot2::xlab("Group") +
    ggplot2::ylab("Number of genes") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = c("#440154FF", "#440154FF", "#2D708EFF", "#2D708EFF")) +
    ggplot2::ggtitle(paste(title, comp, sep = "
"))
  if (is.null(ylim)== FALSE){
    plot <- plot + ggplot2::ylim(ylim)
  }
  plot
  
}



#' GoI plots
#'
#' @param data
#' @param title
#' @param xlab
#'
#' @return plot
GoI_logFC_geom_error_fun <- function(data, title, xlab, ncol) {
  data <- as.data.frame(data ) |>
    dplyr::filter(gene_name %in% target) |>
    dplyr::mutate(FDR_cat = dplyr::case_when(FDR < 0.1 ~ "< 0.1",
                                             FDR >= 0.1 ~"ns",
                                             TRUE  ~"ns")) |>
    dplyr::mutate(gene_name = reorder(gene_name, logFC))
  
  
  
  temp_plot <- ggplot2::ggplot(data) +
    ggplot2::geom_point(ggplot2::aes(x = logFC, y = gene_name, fill = FDR_cat, shape = FDR_cat, size = FDR_cat)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = logFC - lfcSE, xmax = logFC + lfcSE, y = gene_name, width = .18, linetype = FDR_cat)) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::ggtitle(paste(title))+
    ggplot2::facet_wrap(~comp, ncol = ncol)
  
  
  if (dplyr::n_distinct(data$FDR_cat) >1) {
    temp_plot <- temp_plot +
      ggplot2::scale_shape_manual(values = c(23, 5)) +
      ggplot2::scale_size_manual(values = c(3, 2.5)) +
      ggplot2::scale_fill_manual(
        values = c("#440154FF", "#808080"))
  }
  if (dplyr::n_distinct(data$FDR_cat) <=1) {
    temp_plot <-temp_plot +
      ggplot2::scale_shape_manual(values = c(5)) +
      ggplot2::scale_size_manual(values = c(2.5))
  }
  temp_plot
}

#' Data for GoI plots selected by GO pathway
#'
#' @param data
#' @param GO_term
#' @param n
#'
#' @return´
get_top_genes_in_ToI <- function(data, GO_term, n) {
  out <- as.data.frame(data) |>
    dplyr::filter(
      gene_id %in% GO_ToI_fun(GO_term)
    ) |>
    dplyr::slice_min(order_by = FDR, n = n) |>
    dplyr::select(symbol)
  out <- out |> as.vector()
  out <- out[[1]]
  out
}


##------------------------------------------------------------------------
## GSEA pathway stuff
##------------------------------------------------------------------------

#' Plot pathway of interest from GSEA output
#'
#' @param data
#' @param ID
#' @param title
#'
#' @return
plot_POI_fun <- function(data, pathway, title) {
  data <-  data  |>
    dplyr::filter(ID == pathway) |>
    dplyr::mutate(ID= gsub("GOBP_", "", ID),
                  ID= gsub("GOMF_", "", ID),
                  ID = snakecase::to_sentence_case(ID, sep_out=" "),
                  ID = reorder(ID, abs(NES)))
  
  ggplot2::ggplot(data = data, ggplot2::aes(y = ID, x = direction, color = `p.adjust`, size = NES)) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(limits = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(ggplot2::vars(comp), nrow = 1) +
    ggplot2::ylab("") +
    ggplot2::xlab("Direction") +
    ggplot2::ggtitle(title)+
    ggplot2::theme(legend.position = "bottom")
}


#' Count number of DE genes in a specified GO term per group
#'
#' @param data
#' @param GO_term
#'
#' @return
count_DE_genes_in_ToI <- function(data, GO_term) {
  temp <- as.data.frame(data) |>
    dplyr::filter(
      gene_id %in% GO_ToI_fun(GO_term),
      FDR < 0.1
    ) |>
    dplyr::mutate(direction = dplyr::case_when(
      logFC > 0 ~"up",
      logFC < 0 ~ "down"
    )) |> 
    dplyr::group_by(comp, direction, .drop = FALSE) |>
    dplyr::summarise(n = as.numeric(length(gene_id)))
}

#' Bar plot for number of DE genes in a specified GO term per group
#'
#' @param data
#' @param title
#' @param comp
plot_DE_genes_in_ToI <- function(data) {

  data$n <- as.numeric(data$n)
  
  ggplot2::ggplot(data = data) +
    ggplot2::geom_col(ggplot2::aes(x = n, y = comp, fill = direction)) +
    ggplot2::xlab("Group") +
    ggplot2::ylab("Number of genes") +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 2+max(data$n)))+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_manual(values = c("darkred", "midnightblue", "#808080")) +
    ggplot2::geom_text(ggplot2::aes(x = n +2, y = comp, label = n), size = 4, vjust = .5)+
    ggplot2::facet_wrap(~pathway, ncol = 2)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text(size = 12),
                   legend.position = "bottom")
}
