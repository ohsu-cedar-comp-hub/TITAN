#' Topic Elbow Plot
#'
#'
#' This function helps to determine the proper number of topics to use in your LDA model. Plots the rate of perplexity change versus the number of topic. The ideal topic number is the "elbow" of the plot.
#'
#' @param model_dir Directory containing the models created using a varying number of topics.
#' @param SO Seurat object containing the data the model was created with.
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

LDAelbowPlot <- function(model_dir, SO) {
  files <- list.files(path = model_dir, pattern = "Model_")
  # Get model input data
  SO <- FindVariableFeatures(SO, assay = "RNA", nfeatures = 5000)

  # Initialize matrix
  SO.sparse <- GetAssayData(SO, slot = "data", assay = "RNA")

  #Convert to LDA format
  data.use <- SO.sparse[VariableFeatures(SO, assay = "RNA"), ]
  data.use <- data.use * 100
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
