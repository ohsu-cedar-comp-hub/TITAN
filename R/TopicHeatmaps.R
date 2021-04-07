#' Creates vector of annotation colors
#'
#'
#' This function creates a vector of colors similar to the default colors of ggplot2
#'
#' @param n Desired length of vector
#'
#' @examples
#' gg_color_hue(50)
#'
#' @return Vector of colors for plotting
#'
#'
#' @export
#'
#' @import RColorBrewer


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Topic Heatmap
#'
#'
#' This function produces a Heatmap showing how much each cell lies within each topic
#'
#' @param Object Seurat object containing the data the model was created with.
#' @param topics document-topic matrix
#' @param AnnoVector Vector of cell annotations
#' @param AnnoName Name of cell annotation
#' @param clusterTopics Logical whether to cluster topics or not
#'
#' @examples
#' Heatmap(SeuratObj, DocTopMat, SeuratObj$annotation, "NameOfAnnotation")
#'
#' @return pheatmap object
#'
#'
#' @export
#'
#' @import pheatmap
#' @import RColorBrewer

HeatmapTopic <- function(Object,
                    topics,
                    AnnoVector,
                    AnnoName,
                    clusterTopics = F) {

  #Create a dataframe with the annotation information and corresponding colors
  anno_col           <- data.frame(row.names = colnames(Object),
                         Column1=AnnoVector)
  colnames(anno_col) <- AnnoName
  num_colors         <- length(unique(anno_col[,1]))
  anno_colors        <- gg_color_hue(num_colors)
  names(anno_colors) <- sort(unique(anno_col[,1]))
  anno_colors        <- list(Cluster = anno_colors)
  names(anno_colors) <- AnnoName

  #Add annotation color information to topics
  topics <- data.matrix(topics[order(anno_col[,1]),])

  #plot Heatmap
  if (clusterTopics == F) {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  } else {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  }
  return(p1)
}


#' Topic Heatmap - sorted by topic
#'
#'
#' This function produces a Heatmap showing how much each cell lies within each topic. The cells in the Heatmap are sorted based on their scores in a given topic.
#'
#' @param Object Seurat object containing the data the model was created with.
#' @param topics document-topic matrix
#' @param sortByTopic topic to be sorted by
#' @param AnnoVector Vector of cell annotations
#' @param AnnoName Name of cell annotation
#' @param clusterTopics Logical whether to cluster topics or not
#'
#' @examples
#' HeatmapSortByTopic(SeuratObj, DocTopMat, "Topic_16", SeuratObj$annotation, "NameOfAnnotation")
#'
#' @return pheatmap object
#'
#'
#' @export
#'
#' @import pheatmap
#' @import RColorBrewer


HeatmapSortByTopic  <- function(Object,
                                topics,
                                sortByTopic =  "Topic_1",
                                AnnoVector, AnnoName,
                                clusterTopics = F) {

  #Create a dataframe with the annotation information and corresponding colors
  anno_col           <- data.frame(row.names = colnames(Object),
                         Column1=AnnoVector)
  colnames(anno_col) <- AnnoName
  num_colors         <- length(unique(anno_col[,1]))
  anno_colors        <- gg_color_hue(num_colors)
  names(anno_colors) <- sort(unique(anno_col[,1]))
  anno_colors        <- list(Cluster = anno_colors)
  names(anno_colors) <- AnnoName

  #Add annotation color information to topics
  topics <- data.matrix(topics[order(topics[,sortByTopic]),])

  #plot Heatmap
  if (clusterTopics == F) {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1) 
  } else {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1) 
  }
  return(p1)
}


#' Topic Heatmap - sorted by topic and annotation
#'
#'
#' This function produces a Heatmap showing how much each cell lies within each topic. The cells in the Heatmap are sorted based on their scores in a given topic. The annotations are also sorted by their values.
#'
#' @param Object Seurat object containing the data the model was created with.
#' @param topics document-topic matrix
#' @param sortByTopic topic to be sorted by
#' @param AnnoVector Vector of cell annotations
#' @param AnnoName Name of cell annotation
#' @param clusterTopics Logical whether to cluster topics or not
#'
#' @examples
#' HeatmapSortByTopicAsWellAsAnno(SeuratObj, DocTopMat, "Topic_16", SeuratObj$annotation, "NameOfAnnotation")
#'
#' @return pheatmap object
#'
#'
#' @export
#'
#' @import pheatmap
#' @import RColorBrewer


HeatmapSortByTopicAsWellAsAnno  <- function(Object,
                                            topics,
                                            sortByTopic =  "Topic_1",
                                            AnnoVector, AnnoName,
                                            clusterTopics = F) {

  #Create a dataframe with the annotation information and corresponding colors
  anno_col           <- data.frame(row.names = colnames(Object),
                         Column1=AnnoVector)
  colnames(anno_col) <- AnnoName
  num_colors         <- length(unique(anno_col[,1]))
  anno_colors        <- gg_color_hue(num_colors)
  names(anno_colors) <- sort(unique(anno_col[,1]))
  anno_colors        <- list(Cluster = anno_colors)
  names(anno_colors) <- AnnoName

  #Add annotation color information to topics
  topics <- data.matrix(topics[order( anno_col[,1], topics[,sortByTopic]),])

  #plot Heatmap
  if (clusterTopics == F) {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = F,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  } else {
    p1 <- pheatmap(topics,
                   hclustfun = function(x) hclust(x, method="ward.D2"),
                   scale = "row",
                   cluster_cols = T,
                   cluster_rows = F,show_rownames = F,
                   col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
                   annotation_row = anno_col,
                   annotation_names_row = T,
                   annotation_colors = anno_colors,
                   cex=1)
  }
  return(p1)
}



#' #' Clustered Topic Heatmap
#' #'
#' #'
#' #' This function produces a heatmap where the columns are clustered based on the expression pattern of the topics
#' #'
#' #' @param Object Seurat object containing the data the model was created with.
#' #' @param topics document-topic matrix
#' #' @param AnnoVector Vector of cell annotations
#' #' @param AnnoName Name of cell annotation
#' #'
#' #' @examples
#' #' heatmap(SeuratObj, DocTopMat, SeuratObj$annotation, "NameOfAnnotation")
#' #'
#' #' @return pheatmap object
#' #'
#' #'
#' #' @export
#' #'
#' #' @import pheatmap
#' #' @import RColorBrewer
#'
#' heatmapClusterTopics <- function(Object,
#'                                  topics,
#'                                  AnnoVector,
#'                                  AnnoName) {
#'   anno_col <- data.frame(row.names = colnames(Object),
#'                          Column1=AnnoVector)
#'   colnames(anno_col) <- AnnoName
#'   num_colors = length(unique(anno_col[,1]))
#'   anno_colors = gg_color_hue(num_colors)
#'   names(anno_colors) <- sort(unique(anno_col[,1]))
#'   anno_colors <- list(Cluster = anno_colors)
#'   names(anno_colors) <- AnnoName
#'
#'   topics <- data.matrix(topics[order(anno_col[,1]),])
#'
#'   p1 <- pheatmap(topics,
#'                  hclustfun = function(x) hclust(x, method="ward.D2"),
#'                  scale = "row",
#'                  cluster_cols = T,
#'                  cluster_rows = F,show_rownames = F,
#'                  col=colorRampPalette(rev(brewer.pal(11, "RdBu"))[c(1:4,8:11)])(256),
#'                  annotation_row = anno_col,
#'                  annotation_names_row = T,
#'                  annotation_colors = anno_colors,
#'                  cex=1)
#'   return(p1)
#' }
