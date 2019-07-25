
#' Create visual representation of one gene expression according to indexes combinaisons.
#'
#' @import SingleCellExperiment ggplot2
#' @param sce A SingleCellExperiement object.
#' @param goi Gene of interest
#' @return ggplot
#' @export
#' @examples
#' sce <- SingleCellExperiment(assays = list(counts = as.matrix(cleanhop_counts)),
#'                             colData = cleanhop_annotations)
#' indexplot(sce, "B2M")
#' indexplot(sce, "CCR8")

indexplot = function(sce, goi) {
  df <- as.data.frame(colData(sce))
  goi_counts <- NULL
  df$goi_counts <- counts(sce)[goi,]

  plot = ggplot(df,  aes_string( x = "i7_ID"  , y = "i5_ID" )) +
    scale_y_discrete(limits = rev(levels(sce$i5_ID))) +
    geom_point( aes(size = log10(goi_counts+0.1)), shape = 21, stroke = 0.2, color = "darkgrey", fill = "firebrick1", alpha = 0.5) +
    geom_text(aes(label = round(goi_counts,0)))+
    ggtitle(paste0(goi, "_counts"))+
    theme_bw() +
    theme(axis.title.x = element_text(size=12, face="bold"),
          axis.title.y = element_text(size=12, face="bold"),
          axis.text.x = element_text(angle=90),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom" ) +
    facet_grid(Lane_ID ~ Run_ID)

  print(plot)
}
