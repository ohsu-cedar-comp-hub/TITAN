#' Plot2Topics
#'
#'
#' This function allows the user to choose two topics and plot them along the x and y axis of a plot
#'
#' @param SeuratObj Seurat Object containing the data the model was created with.
#' @param model Output of the LDA function
#' @param topics The two topics that will plotted along the axes
#' @param GOI The cells on the plot will be colored by how much they express the given gene. If "CC" is entered, the cells will by colored by their cell cycle phase. Cell cycle must be calculated prior.
#' @param plot_split A column in the metadata of the Seurat Object to split the plots by
#' @param IntIdents Specific factor levels in the plot_split to plot. (If left empty, all factor levels will be plotted)
#'
#' @examples
#' Plot2Topics(SeuratObj = my_SO, model = my_LDA, topics = c(1,2), GOI = "CC")
#'
#' @return Prints plot
#'
#'
#' @export
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#'


plot2Topics <- function(SeuratObj,
                        model,
                        topics,
                        GOI,
                        plot_split,
                        IntIdents) {
  SO <- Seurat_object
  if (length(IntIdents) > 0) {
    keep <- (unname(SO@metadata[,plot_split]) %in% IntIdents)
    SO <- SO[,keep]
  }

  if (plot_split %in% colnames(SO@meta.data)) {
    if (GOI %in% rownames(SO)){
      fp <- FeaturePlot(SO, features = GOI, dims = topics, reduction = "lda", split.by = plot_split, min.cutoff = 'q1', combine = F)
      print(CombinePlots(fp, ncol = 2))
    }
    else {
      if (GOI == "CC"){
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Phase", split.by = plot_split, ncol = 2)
        print(fp)
      }
      else {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, split.by = plot_split, ncol = 2)
        print(fp)
      }
    }
  }
  else {
    if (GOI %in% rownames(SO)){
      fp <- FeaturePlot(SO, features = GOI, dims = topics, reduction = "lda", min.cutoff = 'q1')
      print(fp)
    } else {
      if (GOI == "CC"){
        fp <- DimPlot(SO, reduction = "lda", dims = topics, group.by = "Phase", ncol = 2)
        print(fp)
      }
      else {
        fp <- DimPlot(SO, reduction = "lda", dims = topics, ncol = 2)
        print(fp)
      }
    }
  }
}
