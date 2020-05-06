#' Adds topic information to Seurat metadata
#'
#'
#' This function adds the topic-document information as individual columns to the Seurat metadata for each topic
#'
#' @param model LDA model output
#' @param Object Seurat object containing the data the model was created with.
#'
#' @examples
#' addTopicsToSeuratObject(LDAmodel, SeuratObj)
#'
#' @return Seurat object with updated metadata
#'
#'
#' @export
#'
#' @import Seurat

addTopicsToSeuratObject <- function(model, Object) {
  modelMat           <- t(scale(model$document_expects, center=TRUE, scale=TRUE))
  rownames(modelMat) <- paste(1:ncol(Object), colnames(Object), sep="_")
  colnames(modelMat) <- paste("Topic", 1:ncol(modelMat), sep="_")
  Object@meta.data   <- cbind(Object@meta.data, modelMat)
  return(Object)
} 