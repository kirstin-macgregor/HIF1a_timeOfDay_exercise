##########################################################################
## Project    : Functions for HIF1a analysis
## Data       : RNA-seq
## Author     : Kirstin MacGregor
## Date       : Project start: January 04, 2023. Last update: January 23, 2023
## Version    : v1.0
##########################################################################


##------------------------------------------------------------------------
## Gene ontology
##------------------------------------------------------------------------

#' function for running GSEA and tidying output
#'
#' @param data
#'
#' @return GSEA results
GSEA_fun <- function(data){
  
  comp <- unique(data$comp)
  
  ranked <- data |>
    tidyr::drop_na(entrezid, pvalue, logFC)  |>
    dplyr::mutate(rank = -log10(pvalue) *sign(logFC))  |>
    dplyr::arrange(desc(rank))  |>
    dplyr::pull(rank, entrezid)
  
  term2gene_MF <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF") |>
    dplyr::select(gs_name, entrez_gene)
  term2name_MF <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF") |>
    dplyr::select(gs_name, gs_description) |>
    dplyr::distinct()
  
  GSEA_MF <- clusterProfiler::GSEA(ranked,
                                   TERM2GENE = term2gene_MF,
                                   TERM2NAME = term2name_MF,
                                   pvalueCutoff = 1.00,
                                   minGSSize = 15,
                                   maxGSSize = 500,
                                   verbose=TRUE)
  
  GSEA_MF_out <- tibble::as_tibble(GSEA_MF)  |>
    dplyr::arrange(desc(abs(NES))) |>
    dplyr::mutate(dplyr::across(c("enrichmentScore", "NES"), round, digits=3)) |>
    dplyr::mutate(dplyr::across(c("pvalue", "p.adjust", "qvalue"), scales::scientific),
                  direction= dplyr::case_when(NES>0 ~"Up", NES<0 ~"Down"),
                  comp = paste(comp),
                  pvalue= as.numeric(pvalue),
                  p.adjust= as.numeric(p.adjust),
                  qvalue= as.numeric(qvalue))
  
  term2gene_BP <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP") |>
    dplyr::select(gs_name, entrez_gene)
  term2name_BP <- msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP") |>
    dplyr::select(gs_name, gs_description) |>
    dplyr::distinct()
  
  GSEA_BP <- clusterProfiler::GSEA(ranked,
                                   TERM2GENE = term2gene_BP,
                                   TERM2NAME = term2name_BP,
                                   pvalueCutoff = 1.00,
                                   minGSSize = 15,
                                   maxGSSize = 500,
                                   verbose=TRUE)
  
  GSEA_BP_out <- tibble::as_tibble(GSEA_BP)  |>
    dplyr::arrange(desc(abs(NES))) |>
    dplyr::mutate(dplyr::across(c("enrichmentScore", "NES"), round, digits=3)) |>
    dplyr::mutate(dplyr::across(c("pvalue", "p.adjust", "qvalue"), scales::scientific),
                  direction= dplyr::case_when(NES>0 ~"Up", NES<0 ~"Down"),
                  comp = paste(comp),
                  pvalue= as.numeric(pvalue),
                  p.adjust= as.numeric(p.adjust),
                  qvalue= as.numeric(qvalue))
  
  out <- list("GSEA_MF" = GSEA_MF,
              "GSEA_MF_out" = GSEA_MF_out,
              "GSEA_BP"= GSEA_BP,
              "GSEA_BP_out"= GSEA_BP_out)
}