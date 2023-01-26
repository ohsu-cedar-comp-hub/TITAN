#' Adds imputed Topics to a seruat object
#'
#'
#' This function adds the topic-document information as individual columns to the Seurat metadata for each topic
#'
#' @param model LDA model output
#' @param Object Seurat object containing the data the model was created with.
#' @param TopicPrefix Prefix of Topics to be added to your metaData (Default "ImputedTopic")
#' @param assayName The name of the assay to associate with the reduction
#'
#' @examples
#' ImputeAndAddTopics(SeuratObj, LDAmodel, TopicPrefix = "ImputedTopic")
#'
#' @return Seurat object with updated metadata of imputed geneExpression
#'
#'
#' @export
#'
#' @import Seurat


ImputeAndAddTopics <- function(Object, model, TopicPrefix = "ImputedTopic", assayName = "RNA") {
  Top50Words <- top.topic.words(model$topics, 50, by.score = T)
  wordList <- split(Top50Words, rep(1:ncol(Top50Words), each = nrow(Top50Words)))
  Object <- AddModuleScore(Object, features = wordList, name = paste0(TopicPrefix, "_"))

  Object[["imputedLDA"]] <- CreateDimReducObject(
    embeddings = as.matrix(Object@meta.data %>% select(starts_with(TopicPrefix))),
    key = "imputedLDA_",
    assay = assayName,
    global = TRUE
  )
  return(Object)
}

