#' Topic Elbow Plot
#'
#'
#' This function helps to determine the proper number of topics to use in your LDA model. Plots the rate of perplexity change versus the number of topic. The ideal topic number is the "elbow" of the plot.
#'
#' @param model_dir Directory containing the models created using a varying number of topics.
#' @param Object Object containing the data the model was created with.
#' @param VarFeatures the number of variable features to use in the LDA model. MUST MATCH WITH MODELS IN model_dir
#' @param assayName The name of the assay holding the source data
#'
#' @examples
#' LDAelbowPlot(test_dir, SeuratObj)
#'
#' @return ggplot object
#'
#'
#' @export
#'
#' @import Seurat
#' @import text2vec
#' @import lda
#'
LDAelbowPlot <- function(model_dir,
                         Object, varFeatures = 5000, assayName = "RNA") {
  files <- list.files(path = model_dir, pattern = "Model_")

  # Get model input data
  if (class(Object) == "Seurat") {
    #Normalize and extract the gene expression data from the Seurat Object
    Object        <- NormalizeData(Object, assay = assayName, normalization.method = "CLR")
    Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = varFeatures)
    Object.sparse <- GetAssayData(Object, slot = "data",assay = assayName)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

    #convert data into the proper input format for lda.collapsed.gibbs.sampler
    data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  } else if (class(Object) == "SingleCellExperiment") {
    normalized_sce <- NormalizeData(assay(Object, "counts"), normalization.method = "CLR")
    varFeats <- FindVariableFeatures(normalized_sce)
    varFeats$gene <- rownames(varFeats)
    varFeats <- top_n(varFeats, varFeatures, vst.variance.standardized)

    data.use <- Matrix::Matrix(normalized_sce[varFeats$gene,], sparse = T)
  } else (
    message("Object must be of class singleCellExpriment or Seurat")
  )

  data.use <- data.use * 10
  data.use <- round(data.use)

  #initialize necessary variables
  perp_list     <- NULL
  topic_numbers <- NULL
  RPC           <- NULL
  files         <- files[order(nchar(files), files)]

  for (model_file in files) {
    topic_num     <- as.numeric(gsub("[^0-9]+([0-9]+).*", "\\1", model_file))
    topic_numbers <- c(topic_numbers, topic_num)
    model         <- readRDS(paste0(model_dir, "/", model_file))

    #extract document-term matrix
    docterMat     <- t(as.matrix(data.use))
    docterMat     <- as(docterMat, "sparseMatrix")

    #calculate topic word distribution
    topworddist   <- normalize(model$topics, byrow = T)

    #calculate document topic distribution
    doctopdist    <- normalize(t(model$document_sums), byrow = T)

    #calculate perpelexity
    perp          <- perplexity(docterMat, topworddist, doctopdist)
    perp_list     <- c(perp_list, perp)

    #calculate RPC (rate of perplexity change)
    if (length(perp_list) > 1) {
      RPC_temp <- abs((perp_list[length(perp_list)] - perp_list[length(perp_list) - 1]) / (topic_numbers[length(topic_numbers)] - topic_numbers[length(topic_numbers) - 1]))
      RPC      <- c(RPC, RPC_temp)
    }
  }

  #build plot dataframe and create ggplot object
  plot_df           <- as.data.frame(cbind(topic_numbers[-1], RPC))
  colnames(plot_df) <- c("Topics", "RPC")
  p                 <- ggplot(data = plot_df, aes(x = Topics, y = RPC, group = 1)) + geom_line() + geom_point()

  return(p)
}
