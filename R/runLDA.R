#' Runs LDA Model
#'
#'
#' This function runs an LDA model on scRNA-seq expression data
#'
#' @param Object Object containing the data the model was created with.
#' @param normalizationMethod Normalization method used by Seurat NormalizeData. Options are CLR, LogNormalize and RC.
#' @param ntopics Number of topics to be used in the model. If parallel == TRUE, a vector of topics to run should be inputted
#' @param alpha the value for alpha in the LDA model
#' @param beta the value for beta in the LDA model
#' @param varFeatures the number of variable features to use in the LDA model. The more features that are used, the slower the model will run and the more noise that will be introduced, but the model will be more complete in representing your entire dataset.
#' @param iterations the number of iterations used when learning the LDA model.
#' @param burnin number of iterations to run to allow the model to learn before calculating certain statistics. Models start at random points, so this allows model to get closer to the fit before certain statistics are calculated.
#' @param parallel if TRUE, will run multiple models in parallel. NOT AVAILABLE ON WINDOWS
#' @param outDir if parallel = TRUE, the output directory for the multiple models
#' @param cores Number of cores to use, only applicable if parallel = TRUE
#' @param seed.number random integer to set seed
#' @param assayName The name of the assay holding the source data
#'
#'
#' @examples
#' runLDA(SeuratObj, ntopics = 20)
#'
#' @return LDA Model
#'
#'
#' @export
#'
#' @import Seurat
#' @import lda
#' @import parallel


runLDA <- function(Object,
                   ntopics,
                   alpha = 50,
                   beta = 0.1,
                   varFeatures = 5000,
                   iterations = 500,
                   burnin = 250,
                   seed.number = 8,
                   parallel = F,
                   outDir = NULL,
                   cores = 1,
                   normalizationMethod = "CLR",
                   assayName = "RNA") {

  ## Set seed
  set.seed(seed.number)

  if (class(Object) == "Seurat") {
  #Normalize and extract the gene expression data from the Seurat Object
    Object        <- NormalizeData(Object, assay = assayName, normalization.method = normalizationMethod)
    Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = varFeatures)
    Object.sparse <- GetAssayData(Object, slot = "data",assay = assayName)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
    data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  } else {
    message("Object must be of class Seurat")
  }

  data.use      <- data.use * 10
  data.use      <- round(data.use)
  data.use      <- Matrix::Matrix(data.use, sparse = T)
  sumMat        <- Matrix::summary(data.use)
  cellList      <- split(as.integer(data.use@i),
                    sumMat$j)
  ValueList     <- split(as.integer(sumMat$x),
                     sumMat$j
  )
  cellList      <- mapply(rbind, cellList, ValueList, SIMPLIFY=F)
  Genes         <- rownames(data.use)
  cellList      <- lapply(cellList, function(x) {colnames(x) <- Genes[x[1,]+1];x})

  #Run model
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
    if (parallel) {
      if (!dir.exists(outDir)) {
        dir.create(outDir)
      }
      saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
    } else {
      return(selected.Model)
    }
  }
  if (parallel) {
    mclapply(ntopics, model_maker, mc.cores = cores)
  }
  else {
    Model <- model_maker(ntopics)
    return(Model)
  }
}



