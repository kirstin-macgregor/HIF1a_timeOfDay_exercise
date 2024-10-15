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


## ------------------------------------------------------------------------
## Data exploration
## ------------------------------------------------------------------------

#' Create histogram
#'
#' @param data
#' @param var
#' @param bins
#'
#' @return plot
histo_fun <- function(data, var, bins, title) {
  ggplot2::ggplot(data = data) +
    ggplot2::geom_histogram(ggplot2::aes(x = var), stat = "bin", bins = bins, color = "black", fill = "white") +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("sum counts") +
    ggplot2::theme_bw()
}


#' Summarise count data by sum
#'
#' @param data
#'
#' @return sumamry of counts per column
sum_counts_fun <- function(data) {
  data |>
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE)) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "id") |>
    dplyr::left_join(meta, by = "id")
}

#' pivot longer data to long format
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
dat_long_fun <- function(data) {
  dat_long <- data |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "gene")
  dat_long <- dat_long |>
    tidyr::pivot_longer(!"gene", names_to = "id", values_to = "values") |>
    dplyr::left_join(meta, by = "id")
}

#' Create bar plot showing summed counts per sample
#'
#' @param data
#'
#' @return bar plot
sum_count_plot_fun <- function(data) {
  ggplot2::ggplot(data, ggplot2::aes(y = log2(V1), x = id, fill = genotype)) +
    ggplot2::geom_col(colour = "black") +
    ggplot2::xlab("log2 sum counts") +
    ggplot2::scale_fill_manual(values = c("#482677FF", "#20AF7FFF")) +
    ggplot2::coord_flip()
}

#' Create log counts data frame and status cols for plot
#'
#' @param data
#'
#' @return log counts data and status cols for plot
logCounts_fun <- function(data) {
  logCounts <- log2(data + 1)
  logCounts <- as.data.frame(logCounts)
  logCounts <- logCounts |>
    dplyr::mutate_if(is.numeric, list(~ dplyr::na_if(., Inf))) |>
    dplyr::mutate_if(is.numeric, list(~ dplyr::na_if(., -Inf))) |>
    as.data.frame()
  statusCols <- stringr::str_replace_all(meta$genotype, c(ko = "red", wt = "blue"))
  
  out <- list(logCounts, statusCols)
}

#' Boxplot function
#'
#' @param data datframe
#' @param x x variable
#' @param y y variable
#' @param title plot title for ggtitle
#'
#' @return boxplot
boxplot_function <- function(data, xvar, yvar, title, xlab, ylab) {
  plot <- data |>
    ggplot2::ggplot(ggplot2::aes(x = xvar, y = yvar, fill = genotype)) +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_fill_manual(values = c("#482677FF", "#20AF7FFF")) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' create dataframe with PCA variance explained for each component
#'
#' @param data
#'
#' @return dataframe
var_PCA_fun <- function(data) {
  pca_dat <- stats::prcomp(t(data))
  var_PCA <- summary(pca_dat)
  var_PCA <- as.data.frame(t(var_PCA$importance[2, ]))
  var_PCA <- dplyr::select(var_PCA, 1:10)
}


## ------------------------------------------------------------------------
## PCA
## ------------------------------------------------------------------------

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

#' tidy meta data
#'
#' @param data
#'
#' @return dataframe
meta_dat_tidy_fun <- function(data) {
  data |>
    dplyr::select(c(id, animal_n, set, genotype, batch, zt, exercise, sex, cage, rin)) |>
    dplyr::mutate(
      zt_exercise = paste(zt, exercise, sep = "_"),
      sex_exercise = paste(sex, exercise, sep = "_"),
      genotype_exercise = paste(genotype, exercise, sep = "_"),
      zt_sex = paste(zt, sex, sep = "_"),
      zt_genotype = paste(zt, genotype, sep = "_"),
      genotype_sex = paste(genotype, sex, sep = "_"),
      zt_exercise = forcats::as_factor(zt_exercise),
      sex_exercise = forcats::as_factor(sex_exercise),
      genotype_exercise = forcats::as_factor(genotype_exercise),
      zt_sex = forcats::as_factor(zt_sex),
      zt_genotype = forcats::as_factor(zt_genotype),
      genotype_sex = forcats::as_factor(genotype_sex)
    )
}

#' create dataframe with PCA
#'
#' @param data
#' @param meta
#'
#' @return dataframe
pca_dat_fun <- function(data, meta) {
  pca_dat <- stats::prcomp(t(data))
  pca_dat <- cbind(pca_dat$x[, 1:3], meta[, 1:ncol(meta)]) |>
    as.data.frame() |>
    dplyr::select(!c(animal_n, rin, set)) |>
    tidyr::pivot_longer(!c("PC1", "PC2", "PC3", "id"), names_to = "variable", values_to = "value")
}

#' Create scree plot from PCA data
#'
#' @param data
#'
#' @return plot
scree_plot_fun <- function(data) {
  pca_dat <- stats::prcomp(t(data))
  factoextra::fviz_eig(pca_dat, addlabels = TRUE, ylim = c(0, 50))
}

#' Create PCA plot for PC1 and PC2 with id labels
#'
#' @param data
#' @param variable
#'
#' @return plot
plot_PC1_PC2_label_fun <- function(data, label) {
  temp_dat <- data |> dplyr::filter(variable == "genotype")
  
  temp_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC1, y = PC2, color = value, shape = value)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::xlab(paste("PC1 (", 100 * var_PCA$PC1, " %)", sep = "")) +
    ggplot2::ylab(paste("PC2 (", 100 * var_PCA$PC2, " %)", sep = "")) +
    ggplot2::geom_text(data = temp_dat, ggplot2::aes(x = PC1, y = PC2, label = id)) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste("genotype"))
}

#' Create PCA plot for PC2 and PC3 with labels
#'
#' @param data
#' @param variable
#'
#' @return plot
plot_PC2_PC3_label_fun <- function(data, label) {
  temp_dat <- data |> dplyr::filter(variable == "genotype")
  
  temp_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC2, y = PC3, color = value, shape = value)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2) +
    ggplot2::geom_vline(xintercept = 0, lty = 2) +
    ggplot2::xlab(paste("PC2 (", 100 * var_PCA$PC2, " %)", sep = "")) +
    ggplot2::ylab(paste("PC3 (", 100 * var_PCA$PC3, " %)", sep = "")) +
    ggplot2::geom_text(data = temp_dat, ggplot2::aes(x = PC2, y = PC3, label = id)) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste("genotype"))
}

#' Create PCA plot for PC1 and PC2
#'
#' @param data
#' @param variable
#'
#' @return plot

plot_PC1_PC2_fun <- function(data) {
  out <- list()
  for (i in unique(data$variable)) {
    temp_dat <- data |> dplyr::filter(variable == i)
    
    frame_dat <- temp_dat |>
      dplyr::filter(variable == i) |>
      dplyr::group_by(value) |>
      dplyr::slice(chull(PC1, PC2))
    
    temp_plot <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data = frame_dat, ggplot2::aes(x = PC1, y = PC2, group = value, fill = value, colour = value), alpha = 0.2) +
      if (dplyr::n_distinct(temp_dat$value) < 7) {
        ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC1, y = PC2, color = value, shape = value))
      } else if (dplyr::n_distinct(temp_dat$value) > 7) {
        ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC1, y = PC2, color = value))
      }
    
    temp_plot <- temp_plot + ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::scale_fill_viridis_d() +
      viridis::scale_color_viridis(discrete = TRUE) +
      ggplot2::geom_vline(xintercept = 0, lty = 2) +
      ggplot2::xlab(paste("PC1 (", 100 * var_PCA$PC1, " %)", sep = "")) +
      ggplot2::ylab(paste("PC2 (", 100 * var_PCA$PC2, " %)", sep = "")) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste(i))
    
    out[[paste("PC1_PC2_", i, sep = "")]] <- temp_plot
  }
  plots_out <- out
}

#' Create PCA plot for PC2 and PC3
#'
#' @param data
#' @param variable
#'
#' @return plot
plot_PC2_PC3_fun <- function(data) {
  out <- list()
  for (i in unique(data$variable)) {
    temp_dat <- data |> dplyr::filter(variable == i)
    
    frame_dat <- temp_dat |>
      dplyr::filter(variable == i) |>
      dplyr::group_by(value) |>
      dplyr::slice(chull(PC2, PC3))
    
    temp_plot <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data = frame_dat, ggplot2::aes(x = PC2, y = PC3, group = value, fill = value, colour = value), alpha = 0.2) +
      if (dplyr::n_distinct(temp_dat$value) < 7) {
        ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC2, y = PC3, color = value, shape = value))
      } else if (dplyr::n_distinct(temp_dat$value) > 7) {
        ggplot2::geom_point(temp_dat, mapping = ggplot2::aes(x = PC2, y = PC3, color = value))
      }
    temp_plot <- temp_plot + ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::scale_fill_viridis_d() +
      viridis::scale_color_viridis(discrete = TRUE) +
      ggplot2::geom_vline(xintercept = 0, lty = 2) +
      ggplot2::xlab(paste("PC2 (", 100 * var_PCA$PC2, " %)", sep = "")) +
      ggplot2::ylab(paste("PC3 (", 100 * var_PCA$PC3, " %)", sep = "")) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste(i))
    
    out[[paste("PC2_PC3_", i, sep = "")]] <- temp_plot
  }
  plots_out <- out
}