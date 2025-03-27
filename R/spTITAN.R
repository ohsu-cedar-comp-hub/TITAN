#' Finds centroids 
#'
#'
#' This function finds centroids in the x,y space and plots them to visualize their location
#' 
#' @param Object Object containing the data the model was created with.
#' @param pos_df Dataframe containing spatial x and y coordinates of the cells. The rownames of the dataframe need to be the cell names, followed by the x coordinates in column1, then the y coordinates in column2.
#' @param k Number of centroids to use
#' @param seed.number random integer to set seed
#'
#'
#' @examples
#' getCentroids(SeuratObj, cell_positions, k)
#'
#' @return dataframe
#'
#'
#' @export
#'
#' @import Seurat
#' @import stats
#' @import ggplot2

getCentroids <- function(Object,
                         pos_df,
                         k,
                         seed.number = 8) {
    
    set.seed(seed.number)

    pos_df <- pos_df[colnames(Object), 1:2]
    pos_df$points <- "data"
    
    km <- kmeans(pos_df[,1:2], k, nstart = 25)

    center_df <- as.data.frame(km$center)
    center_df$points <- "centroid"
    rownames(center_df) <- paste0("centroid_", 1:nrow(center_df))

    total_df <- rbind(pos_df, center_df)
    colnames(total_df) <- c("x_pos", "y_pos", "points")

    p <- ggplot(total_df, mapping = aes(x=x_pos, y=y_pos, color = points)) + geom_point() + theme_classic()
    print(p)
    return(total_df)
}



#' Calculates input matrix 
#'
#'
#' This function calculates the scaled inverse distance matrix and concatenates it with the gene matrix in preparation for inputting it into the spTITAN algorithm
#' 
#' @param Object Seurat object containing gene expression information.
#' @param centroid_df Dataframe output by getCentroids. Inverse distance matrix is calculated on this matrix
#' @param quant_cut Quantile value for the number of centroid distance values to keep for each cell. Recommended range is between .75 to .95, but other values are viable.
#' @param scale_factor Integer value to scale down distance values by. Recommended range is between 3 to 15, but otherr values are viable
#' @param seed.number random integer to set seed
#'
#'
#' @examples
#' calcDistMat(Object, centroid_df)
#'
#' @return dataframe
#'
#'
#' @export
#'
#' @import Seurat
#' @import stats
#' @import ggplot2
#' @import pdist


calcDistMat <- function(Object,
                        centroid_df,
                        quant_cut,
                        scale_factor,
                        seed.number = 8) {

    set.seed(seed.number)

    pos_df <- subset(centroid_df, points == "data")
    center_df <- subset(centroid_df, points == "centroid")

    dist_df <- as.matrix(pdist(pos_df[,1:2], center_df[,1:2]))
    rownames(dist_df) <- rownames(pos_df)
    colnames(dist_df) <- rownames(center_df)

    inv_dists <- max(dist_df) - dist_df

    test <- t(inv_dists)

    for (i in 1:ncol(test)) {
        temp <- test[,i]
        cutoff <- quantile(temp, probs = quant_cut)
        temp[which(temp < cutoff)] <- 0
        test[,i] <- temp
    }

    test <- test/scale_factor

    Object.sparse <- GetAssayData(SO, slot = "counts", assay = "RNA")
    final_mat <- rbind(Object.sparse, test)

    return(final_mat)
}


#' Runs LDA Model
#'
#'
#' This function runs an LDA model on scRNA-seq expression data
#'
#' @param input_df Output of calcDistMat that is being fed into spTITAN algorithm
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
#' spTITAN(input_df, ntopics = 20)
#'
#' @return LDA Model
#'
#'
#' @export
#'
#' @import Seurat
#' @import lda
#' @import parallel
#' @import stats


spTITAN <- function(input_df,
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
                    assayName = "RNA",
                    visualizeCentroids = TRUE) {

  ## Set seed
  set.seed(seed.number)


    temp_SO <- CreateSeuratObject(counts=input_df)
    Object        <- NormalizeData(temp_SO, assay = assayName, normalization.method = normalizationMethod)
    Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = varFeatures)
    
    Object.sparse <- GetAssayData(Object, slot = "data",assay = assayName)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
    data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  
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


#' spTITAN Topic Elbow Plot
#'
#'
#' This function helps to determine the proper number of topics to use in your LDA model. Plots the rate of perplexity change versus the number of topic. The ideal topic number is the "elbow" of the plot.
#'
#' @param model_dir Directory containing the models created using a varying number of topics.
#' @param input_df Output of calcDistMat
#' @param VarFeatures the number of variable features to use in the LDA model. MUST MATCH WITH MODELS IN model_dir
#' @param assayName The name of the assay holding the source data
#'
#' @examples
#' spTITANelbowPlot(test_dir, SeuratObj)
#'
#' @return ggplot object
#'
#'
#' @export
#'
#' @import Seurat
#' @import text2vec
#' @import lda



spTITANelbowPlot <- function(model_dir,
                         input_df, varFeatures = 5000, assayName = "RNA") {
  files <- list.files(path = model_dir, pattern = "Model_")

  # Get model input data
  temp_SO <- CreateSeuratObject(counts=input_df)
    #Normalize and extract the gene expression data from the Seurat Object
    Object        <- NormalizeData(temp_SO, assay = assayName, normalization.method = "CLR")
    Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = varFeatures)
    
    Object.sparse <- GetAssayData(Object, slot = "data",assay = assayName)
    Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

    #convert data into the proper input format for lda.collapsed.gibbs.sampler
    data.use      <- Matrix::Matrix(Object.sparse, sparse = T)
  
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
