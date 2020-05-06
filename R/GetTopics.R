#' Extract Topic Information
#'
#'
#' This function extracts the topic information from the LDA model. The output matrix is the document-topic matrix showing how much each cell corresponds with each topic
#'
#' @param model LDA model.
#' @param Object Seurat object containing the data the model was created with.
#'
#' @examples
#' GetTopics(LDA_model, SeuratObj)
#'
#' @return document-topic matrix
#'
#'
#' @export


GetTopics <- function(model, Object) {
  modelMat           <- t(scale(model$document_expects, center=TRUE, scale=TRUE))
  rownames(modelMat) <- colnames(Object)
  colnames(modelMat) <- paste("Topic", 1:ncol(modelMat), sep="_")
  return(modelMat)
}