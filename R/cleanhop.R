
#' Remove potential misassigned reads due to index hopping in single cell RNAseq data generated with combinatorial dual indexes.
#'
#' @import SingleCellExperiment
#' @param sce A SingleCellExperiement object.
#' @param frct_value Fraction of total reads to substract. Default value = 0.005.
#' @return \code{sce} A SingleCellExperiement object.
#' @export
#' @examples
#' sce <- SingleCellExperiment(assays = list(counts = as.matrix(cleanhop_counts)),
#'                             colData = cleanhop_annotations)
#' cleaned_sce = cleanhop(sce)
#' cleaned_sce = cleanhop(sce, frct_value = 0.01)

cleanhop = function(sce, frct_value = 0.005) {
  for (this_run in levels(sce$Run_ID)){
    for (this_lane in levels(sce[,sce$Run_ID == this_run]$Lane_ID)) {
      for (this_i5 in levels(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]$i5_ID)) {
        cells_this_i5 <- sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]$i5_ID == this_i5
        this_i5_gene_count_sum <- rowSums(counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i5]))
        values_to_substract <- round(this_i5_gene_count_sum * frct_value)
        counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i5]) <- counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i5]) - values_to_substract
      }
      counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane])[counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]) < 0] <- 0

      for (this_i7 in levels(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]$i7_ID)) {
        cells_this_i7 <- sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]$i7_ID == this_i7
        this_i7_gene_count_sum <- rowSums(counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i7]))
        values_to_substract <- round(this_i7_gene_count_sum * frct_value)
        counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i7]) <- counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane][, cells_this_i7]) - values_to_substract
      }
      counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane])[counts(sce[,sce$Run_ID == this_run & sce$Lane_ID == this_lane]) < 0] <- 0
    }
  }

  return(sce)
}



#' Example single cell RNAseq data reduced to 8 genes and 768 cells, counts
#'
#' A dataset containing single cell RNAseq gene counts for 8 genes in 768 cells barcoded with Nextera XT indexes
#'
#' @format a data frame with 8 row and 768 columns:
#' \describe{
#'     \item{columns}{Cell_ID}
#' }
"cleanhop_counts"

#' Example single cell RNAseq data reduced to 8 genes and 768 cells, column annotations
#'
#' A dataset containing single cell RNAseq gene annotation for 768 cells barcoded with Nextera XT indexes
#'
#' @format a data frame with 768 row and 4 columns:
#' \describe{
#'     \item{Run_ID}{Sequencing run associated with expression data (counts) for each cell}
#'     \item{Lane_ID}{Sequencing lane associated with expression data (counts) for each cell}
#'     \item{i5_ID}{First index associated to each cell}
#'     \item{i7_ID}{Second index associated to each cell}
#'     \item{row}{Cell_ID}
#' }
"cleanhop_annotations"

#' Counts table generate after the use of cleanhop function on example single cell RNAseq data reduced to 8 genes and 768 cells
#'
#' A dataset containing single cell RNAseq gene counts for 8 genes in 768 cells barcoded with Nextera XT indexes
#'
#' @format a data frame with 8 row and 768 columns:
#' \describe{
#'     \item{columns}{Cell_ID}
#' }
"cleaned_cleanhop_counts"




# sce = readRDS("/Users/OBricard/Documents/R/cleanhop/data/sce.rds")
# cleanhop_counts = as.data.frame(counts(sce))
# usethis::use_data(cleanhop_counts)
# cleanhop_annotations = as.data.frame(colData(sce))
# usethis::use_data(cleanhop_annotations)
# cleaned_sce = cleanhop(sce)
# cleaned_cleanhop_counts = as.data.frame(counts(cleaned_sce))
# usethis::use_data(cleaned_cleanhop_counts)

