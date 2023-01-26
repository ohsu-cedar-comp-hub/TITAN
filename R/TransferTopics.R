#' Transfers topics from one dataset to another
#'
#'
#' This function takes the topics from one dataset and, using the word-topic distribution from that dataset, to impute these topics onto a new dataset.
#'
#' @param Object Seurat object to transfer topics onto
#' @param Model LDA model using data from another dataset
#' @param TopicPrefix Common prefix for topics in the LDA model
#' @param assayName The name of the assay to associate with the reduction
#' 
#'
#' @examples
#' TransferTopics(NewSeuratObj, Old_model)
#'
#' @return Seurat object with topics transferred into its metadata
#'
#' @export
#'
#' @import Seurat


TransferTopics <- function(Object, Model, TopicPrefix = "Topic", assayName = "RNA") {
  #extracts the top genes for each topic from the word-topic distribution
  Top50Words <- top.topic.words(Model$topics, 50, by.score = T)
  wordList   <- split(Top50Words, rep(1:ncol(Top50Words), each = nrow(Top50Words)))
  
  #Creates a module score in the new dataset using those top genes
  SeuratObject <- AddModuleScore(SeuratObject, features = wordList, name = paste0(TopicPrefix, "_"))
  
  SeuratObject[["imputedLDA"]] <- CreateDimReducObject(
    embeddings = as.matrix(SeuratObject@meta.data %>% select(starts_with(TopicPrefix))),
    key = "imputedLDA_",
    assay = assayName,
    global = TRUE
  )
  return(SeuratObject)
}