#' TopicWordTable function for SingleCellExperiment
#'
#'
#' Dataframe with the top n genes per topic
#'
#' @param Model  LDA Model Object
#' @param nGenes Number of genes do you want returned per topic
#'
#' @examples
#' SCE_FeaturePlot(Model = LDAmodel, nGenes = 50)
#'
#' @return Returns a dataframe with the top n genes per topic
#'
#'
#' @export
#'
#' @import lda


TopicWordTable <- function(Model, nGenes = 50) {
  LDA_Topics <- top.topic.words(Model$topics, num.words = nGenes, by.score = T)
  colnames(LDA_Topics) <- paste("Topic", 1:ncol(LDA_Topics), sep = "_")
  return(LDA_Topics)
}


