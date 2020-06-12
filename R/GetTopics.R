#' Extract Topic Information
#'
#'
#' This function extracts the cell-topic information from the LDA model. The output matrix is the cell-topic matrix giving a normalized score of how much each cell corresponds with each topic
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


GetTopics <- function(model,
                      Object) {
  ## Get Cell Topic Scores and scale across Topics
  modelMat           <- t(scale(model$document_expects, center=TRUE, scale=TRUE))
  rownames(modelMat) <- colnames(Object)
  colnames(modelMat) <- paste("Topic", 1:ncol(modelMat), sep="_")

  ## Return cell-document matrix
  return(modelMat)
}


#' Extract Top Topic Genes
#'
#'
#' This function extracts the top n scoring genes for each gene
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

  ## LDA function for extracting top n genes
  Topic_Genes        <- lda::top.topic.words(model$topics, num.words = ngenes, by.score = T)
  colnames(Topic_Genes) <- paste("Topic", 1:ncol(Topic_Genes), sep = "_")

  ##Output as a table
  return(Topic_Genes)
}
