#' Runs a series of LDA models in parallel
#'
#'
#' This function runs a series of LDA models on scRNA-seq expression data in parallel.
#'
#' @param Object Seurat object containing the data the model was created with.
#' @param ntopics Number of topics to be used in the model
#' @param outdir Directory for all of the output models
#' @param alpha the value for alpha in the LDA model
#' @param beta the value for beta in the LDA model
#' @param varFeatures the number of variable features to use in the LDA model. The more features that are used, the slower the model will run and the more noise that will be introduced, but the model will be more complete in representing your entire dataset.
#' @param iterations the number of iterations used when learning the LDA model.
#' @param burnin number of iterations to run to allow the model to learn before calculating certain statistics. Models start at random points, so this allows model to get closer to the fit before certain statistics are calculated.
#' @param top_start Lowest number of topics to run a model with
#' @param top_end Highest number of topics to run a model with
#'
#'
#' @examples
#' parallel_runLDA(SeuratObj, ntopics = 20, "model_dir/")
#'
#' @return saves the different models to the given output directory
#'
#'
#' @export
#'
#' @import Seurat
#' @import lda
#' @import parallel


parallel_runLDA <- function(Object, ntopics, outdir, alpha = 50, beta = 0.1, varFeatures = 5000, iterations = 500, burnin = 250, top_start = 5, top_end = 50) {

  #Normalize and extract the gene expression data from the Seurat Object
  Object        <- NormalizeData(Object, assay = "RNA", normalization.method = "CLR")
  Object        <- FindVariableFeatures(Object, assay = "RNA", nfeatures = varFeatures)
  Object.sparse <- GetAssayData(Object[VariableFeatures(Object, assay = "RNA"),], slot = "data",assay = "RNA")
  #Object.sparse <- Object.sparse[VariableFeatures(Object, assay = "RNA"),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
  data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  data.use      <- round(data.use * 10)
  #data.use      <- Matrix::Matrix(data.use, sparse = T)
  sumMat        <- Matrix::summary(data.use)
  cellList      <- split(as.integer(data.use@i),
                         sumMat$j)
  ValueList     <- split(as.integer(sumMat$x),
                         sumMat$j
  )
  cellList      <- mapply(rbind, cellList, ValueList, SIMPLIFY=F)
  Genes         <- rownames(data.use)
  cellList      <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1,]+1];x})

  #create vector of topic numbers
  topic_options <- seq(top_start, top_end, 5)

  #Run models using the different numbers of topics in parallel
  model_maker <- function(topics) {
    selected.Model <- lda.collapsed.gibbs.sampler(
      cellList,
      topics,
      Genes,
      num.iterations=iterations,
      alpha=alpha,
      eta=beta,
      compute.log.likelihood=TRUE,
      burnin=burnin)[-1]
    saveRDS(selected.Model, paste0(outdir, "Model_", as.character(topics), "topics.rds"))
  }

  mclapply(topic_options, model_maker, mc.cores = length(topic_options))
}



