#' Adds topic information to Seurat metadata
#'
#'
#' This function adds the topic-document information as individual columns to the Seurat metadata for each topic
#'
#' @param model LDA model output
#' @param Object Seurat object containing the data the model was created with.
#' @param assayName The name of the assay to associate with the reduction
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

addTopicsToSeuratObject <- function(model,
                                    Object, assayName = "RNA") {

  ## Get Cell Topic Scores and scale across Topics
  modelMat           <- t(scale(model$document_expects, center=TRUE, scale=TRUE))
  rownames(modelMat) <- paste(1:ncol(Object), colnames(Object), sep="_")
  colnames(modelMat) <- paste("Topic", 1:ncol(modelMat), sep="_")

  ## Add to metaData
  Object@meta.data   <- cbind(Object@meta.data, modelMat)

  ## Add lda topics to a dim reduc
  rownames(modelMat) = colnames(Object)
  Object[["lda"]] <- CreateDimReducObject(
    embeddings = modelMat,
    key = "lda_",
    assay = assayName,
    global = TRUE
  )
  return(Object)
}


