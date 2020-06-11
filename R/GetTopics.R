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


#' Extract Top Topic Genes
#'
#'
#' This function extracts the top scoring genes for each gene
#'
#' @param model LDA model.
#' @param ngenes Number of genes to extract for each topic.
#' 
#'
#' @examples
#' TopTopicGenes(LDA_model, 50)
#'
#' @return matrix containing top genes for each topic
#'
#'
#' @export
#' 
#' @import lda


TopTopicGenes <- function(model, ngenes) {
  Topic_Genes        <- lda::top.topic.words(model$topics, num.words = ngenes, by.score = T)
  colnames(Topic_Genes) <- paste("Topic", 1:ncol(Topic_Genes), sep = "_")
  return(Topic_Genes)
}