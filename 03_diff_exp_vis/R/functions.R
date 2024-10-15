##########################################################################
## Project    : Functions for HIF1a analysis
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


## ------------------------------------------------------------------------
## Visualisation
## ------------------------------------------------------------------------

#' Volcano plot
#'
#' @param data
#'
#' @return a volcano plot
volcano_fun <- function(data) {

  ggplot2::ggplot(data = data, ggplot2::aes(x = logFC, y = -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(colour = dir_cat),
      size = 1, alpha = 1.0) +
    ggplot2::scale_colour_manual(values = c("darkred", "darkblue", "darkgrey")) +
    ggrepel::geom_text_repel(data = data |>  dplyr::group_by(comp) |> dplyr::filter(FDR < 0.1) |> dplyr::slice_min(FDR, n= 5, with_ties = FALSE), 
                             ggplot2::aes(label = gene_name), box.padding = 0.5, max.overlaps = Inf, size = 3.5) +
    ggplot2::geom_vline(xintercept = c(-1, 0, 1), linetype = "dotted") +
    ggplot2::theme_bw() +
    ggplot2::xlim(values=c(-10,5))+
    ggplot2::facet_wrap(~comp)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   strip.text = ggplot2::element_text(size = 12))
}

#' tidy data for upset plot
#'
#' @param data
#'
#' @return data arranged for Mladens upset plot function
upset_plot_dat_fun <- function(data) {
  
  dat_temp <- data |> 
    dplyr::group_by(comp) |> 
    dplyr::group_split() |> 
    setNames(sort(unique(data$comp)))
  
  out <- lapply(dat_temp, function(data){
    data |> 
      dplyr::filter(FDR < 0.1) |>
    dplyr::select(gene_id) |>
    dplyr::rename(comp_name = gene_id) |>
    t() |>
    as.vector() 
  }
  )
  out
}

#' Make list of data with columns needed for correlation plot
#'
#' @param data
#'
#' @return plot
make_corr_dat_fun <- function(data) {
  corr_dat <- data |> dplyr::select(logFC, FDR, gene_name)
  corr_dat
}

#######################################
###              Theme              ###
#######################################
library(ggplot2)
library(ggsci)

myTheme <- list(
  "bw" = theme(
    panel.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "black", size = .5),
    axis.ticks = ggplot2::element_line(colour = "black", size = .5),
    text = ggplot2::element_text(size = 8, colour = "black", family = "sans"),
    plot.title = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    axis.text = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    axis.title = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    strip.text = ggplot2::element_text(size = 12, colour = "black", family = "sans", margin = ggplot2::margin(2,2,2,2,"pt")),
    legend.text = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    legend.title = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    legend.key.height = ggplot2::unit(.75, "line"),
    legend.key.width = ggplot2::unit(.75, "line"),
    legend.background = ggplot2::element_blank()
  ),
  
  "min" = theme(
    panel.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 8, colour = "black", family = "sans"),
    plot.title = ggplot2::element_text(size = 9, colour = "black", family = "sans"),
    axis.text = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    axis.title = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    strip.text = ggplot2::element_text(size = 12, colour = "black", family = "sans", margin = ggplot2::margin(2,2,2,2,"pt")),
    legend.text = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    legend.title = ggplot2::element_text(size = 12, colour = "black", family = "sans"),
    legend.key.height = ggplot2::unit(.75, "line"),
    legend.key.width = ggplot2::unit(.75, "line"),
    legend.background = ggplot2::element_blank()
  )
)

### Palette ###
col <- ggsci::pal_npg()((10))

#######################################
###          Upset from List        ###
#######################################
upset_gg <- function(setList, setOrder = NULL, intersectOrder = "level", color = "black"){
  ### setList - named list of sets to intersect, sets in list are vectors of elements
  ### setOrder - order of sets on graph rows, default NULL sets to order of appearence
  ### intersectOrder - order of intersects in columns
  ###                  (default) level returns in order of appearance and combinations
  ###                  freq returns by largest number of intecepts
  ### color - color of bars and dots in the plot, default black
  
  library(ggplot2)
  library(patchwork)
  
  ################
  ### Set Plot ###
  ################
  ### Plot show total number of hits per set ###
  ### Plotting length of each element in list ###
  sPlot <- data.frame(
    "names" = names(sapply(setList, length)),
    "size" = sapply(setList, length)
  )
  
  if(!is.null(setOrder)){
    sPlot$names <- factor(sPlot$names, levels = setOrder)
  } else if (is.null(setOrder)){
    sPlot$names <- factor(sPlot$names, levels = names(setList))
  }
  
  sPlot <- droplevels(sPlot[sPlot$size != 0,])
  
  sP <- ggplot2::ggplot(sPlot, ggplot2::aes(x = size, y = names))+
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = names))+
   # ggplot2::geom_bar(stat = "identity", fill = color)+
    ggplot2::scale_fill_manual(values= c("darkgreen", "lightgreen", "darkgreen", "lightgreen")) +
    ggplot2::xlab("Set Size (n)")+
    ggplot2::scale_x_reverse()+
    ggplot2::theme_minimal()+
    myTheme$min+
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = element_text(angle = 45, hjust=1),
                   legend.position = "none")
  
  
  #########################
  ### Intersection Plot ###
  #########################
  # Combination of all setList sets (i.e. groups to intersect)
  # Combinations in reverse order, so intersects of all groups take priority of unique changes
  comb <- list()
  
  for(i in 1:length(setList)){
    comb[[i]] <- combn(names(setList), i, simplify = FALSE)
  }
  
  comb <- rev(unlist(comb, recursive = FALSE))
  
  # Create Set Intersects
  # Dropping duplicates detected in previous Intersects
  ints <- list()
  
  for(i in 1:length(comb)){
    ints[[i]] <- Reduce(intersect, setList[comb[[i]]])
    
    if(i > 1){
      ints[[i]] <- setdiff(ints[[i]], unlist(ints[1:(i-1)]))
    }
  }
  
  # Create Data Frame
  # Intersect ordering set by intersectOrder - either level or by frequency
  iPlot <- data.frame(
    "intersect" = paste("int", 1:length(comb), sep = ""),
    "sets" = unlist(lapply(comb, function(x) paste(x, collapse = "#sep#"))),
    "size" = lengths(ints)
  )
  
  if(intersectOrder == "level"){
    iPlot$sets <- factor(iPlot$sets, levels = rev(iPlot$sets))
  } else if (intersectOrder == "freq"){
    iPlot$sets <- factor(iPlot$sets, levels = iPlot$sets[order(iPlot$size, decreasing = TRUE)])
  }
  
  iPlot$intersect <- factor(iPlot$intersect, levels = iPlot$intersect[order(iPlot$sets)])
  
  iPlot <- iPlot[iPlot$size != 0,]
  
  iP <- ggplot2::ggplot(iPlot, ggplot2::aes(intersect, size))+
    #ggplot2::geom_bar(stat = "identity", fill = color)+
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill=intersect))+
    ggplot2::scale_fill_manual(values= c("darkgreen", "lightgreen", "darkgreen", "lightgreen", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"))+
    ggplot2::geom_text(ggplot2::aes(y = size * 1.1, label = size), size = 4, vjust = 0, color = color)+
    ggplot2::ylab("Intersect size (n)")+
    ggplot2::theme_minimal()+
    myTheme$min+
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "none")
  
  ###################
  ### Matrix Plot ###
  ###################
  # Show which sets are connected by intersections
  # Reordering levels to correspond to the set plot and intersect plot
  mPlot <- strsplit(as.character(iPlot$sets), split = "#sep#")
  mPlot <- data.frame(
    "intersect" = rep(as.character(iPlot$intersect), lengths(mPlot)),
    "sets" = unlist(mPlot)
  )
  
  # Add grey dots
  mPlot$group <- mPlot$intersect
  mPlot$col <- "a"
  
  grey <- expand.grid("intersect" = unique(mPlot$intersect), "sets" = unique(mPlot$sets))
  grey <- grey[!(interaction(grey$intersect, grey$sets) %in% interaction(mPlot$intersect, mPlot$sets)), ]
  grey$group <- 1:nrow(grey)
  grey$col <- "b"
  
  mPlot <- rbind(mPlot, grey)
  
  # Order corresponding to other plots
  mPlot$intersect <- factor(mPlot$intersect, levels = levels(iPlot$intersect))
  mPlot$sets <- factor(mPlot$sets, levels = levels(sPlot$names))
  
  # Adding row highlight (different for odd and even rows)
  mPlot$row <- "a"
  mPlot$row[as.numeric(mPlot$sets)%%2 != 0] <- "b"
  
  mP <- ggplot2::ggplot(mPlot, ggplot2::aes(intersect, sets, group = group, color = col))+
    ggplot2::geom_tile(aes(fill = row), color = "#00000000", alpha = .2, show.legend = FALSE)+
    ggplot2::geom_point(size = 3, show.legend = FALSE)+
    ggplot2::geom_line(size = 1, show.legend = FALSE)+
    ggplot2::scale_color_manual(values = c(color, "grey"))+
    ggplot2::scale_fill_manual(values = c("#00000000", "grey80"))+
    ggplot2::theme_minimal()+
    myTheme$min+
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())
  
  ################################
  ### Assemble Plot and Output ###
  ################################
  # Assemble plot using patchwork
  lCol <- ggplot2::theme(plot.margin = ggplot2::margin(r = 0))
  rCol <- ggplot2::theme(plot.margin = ggplot2::margin(l = 0))
  
  upset <- invisible(plot_spacer() + lCol + iP + rCol + sP + lCol + mP + rCol +
                       plot_layout(ncol = 2, byrow = TRUE, widths = c(.2,.8), heights = c(.8,.2)))
  
  outList <- invisible(
    list(
      "upset" = upset,
      "intersect_plot" = iP,
      "set_plot" = sP,
      "matrix_plot" = mP,
      "data_frames" = list(
        "intersect" = iPlot,
        "set" = sPlot,
        "matrix" = mPlot
      )
    )
  )
  
  return(outList)
}