#' Topic Elbow Plot
#'
#'
#' This function helps to determine the proper number of topics to use in your LDA model. Plots the rate of perplexity change versus the number of topic. The ideal topic number is the "elbow" of the plot.
#'
#' @param model_dir Directory containing the models created using a varying number of topics.
#' @param SO Object containing the data the model was created with.
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
#' @import LICORS
#' @import lda
#' @import SingleCellExperiment

LDAelbowPlot <- function(model_dir,
                         SO) {
  files <- list.files(path = model_dir, pattern = "Model_")

  # Get model input data
  if (class(Object) == "Seurat") {
    #Normalize and extract the gene expression data from the Seurat Object
    Object        <- NormalizeData(Object, assay = "RNA", normalization.method = "CLR")
    Object        <- FindVariableFeatures(Object, assay = "RNA", nfeatures = varFeatures)
    Object.sparse <- GetAssayData(Object, slot = "data",assay = "RNA")
    Object.sparse <- Object[VariableFeatures(Object, assay = "RNA"),]

    #convert data into the proper input format for lda.collapsed.gibbs.sampler
    data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  } else if (class(Object) == "SingleCellExperiment") {
    normalized_sce <- NormalizeData(assay(Object, "counts"), normalization.method = "CLR")
    varFeats <- FindVariableFeatures(normalized_sce)
    varFeats$gene <- rownames(varFeats)
    varFeats <- top_n(varFeats, 5000, vst.variance.standardized)

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
    topworddist   <- LICORS::normalize(model$topics, byrow = T)

    #calculate document topic distribution
    doctopdist    <- LICORS::normalize(t(model$document_sums), byrow = T)

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
