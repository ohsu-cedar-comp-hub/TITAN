#' FeaturePlot function for SingleCellExperiment
#'
#'
#' Plots the Seurat FeaturePlot function for the topics within SingleCellExperiment objects
#'
#' @param Object SingleCellExperiment object with topics added to it
#' @param pt.size Size of points on plot
#' @param topics The topics to plot on the feature plot
#' @param ncol Number of columns to use when faceting plots
#' @param min.cutoff Vector of minimum and maximum cutoff values for each feature
#'
#' @examples
#' SCE_FeaturePlot(SCE, 0.01, 1:12, 4, 'q5')
#'
#' @return Returns a ggplot object or list of ggplot objects
#'
#'
#' @export
#'
#' @import Seurat
 

SCE_FeaturePlot <- function(object, pt.size, topics, ncol, min.cutoff) {
  temp_SO <- as.Seurat(object, counts = names(T47D_SCE@assays)[1])
  FeaturePlot(temp_SO, pt.size = pt.size, features = paste("Topic", topics, sep = "_"), ncol = ncol, min.cutoff = min.cutoff)
}


