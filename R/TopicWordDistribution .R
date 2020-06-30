#' TopicWordTable function for SingleCellExperiment
#'
#'
#' Returns a dataframe with the top n genes per topic
#'
#' @param Model  LDA Model Object
#' @param nGenes Number of genes do you want returned per topic
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


SCE_FeaturePlot <- function(Model, nGenes = 50) {
  LDA_Topics <- top.topic.words(Model$topics, num.words = nGenes, by.score = T)
  colnames(LDA_Topics) <- paste("Topic", 1:ncol(LDA_Topics), sep = "_")
  return(LDA_Topics)
}


