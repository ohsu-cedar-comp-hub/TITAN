args <- commandArgs()

help <- function(){
    cat("compare_corr.R :
Compares accessibility/RNA correlation and methylation/RNA correlation\n")
    cat("Usage: \n")
    cat("--rna         : RNA object (Seurat)                                       [required]\n")
    cat("--spatLocs    : matrix containing spatial cell locations (rds)            [required]\n")
    cat("--ntopics     : number of topics to use                                   [required]\n")
    cat("--distMat     : pre-calculated inverse distance matrix                    [required]\n")
    cat("--outdir      : output directory for models                               [required]\n")
    cat("\n")
    q()
}

io   <- list()
opts <- list()

## Save values of each argument
if( !is.na(charmatch("--help",args)) || !is.na(charmatch("--help",args)) ){
    help()
} else {
    io$rna         <- sub( '--rna=', '', args[grep('--rna', args)] )
    io$spatLocs    <- sub( '--spatLocs=', '', args[grep('--spatLocs', args)] )
    io$outdir      <- sub( '--outdir=', '', args[grep('--outdir', args)] )
    io$distMat     <- sub( '--distMat=', '', args[grep('--distMat', args)] )
    opts$ntopics   <- sub( '--ntopics=', '', args[grep('--ntopics', args)] )
}

library(Seurat)
library(TITAN)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(parallel)
library(lda)
library(text2vec)
library(LICORS)
library(reticulate)

scipy_sparse <- import("scipy.sparse")

csr_matrix <- scipy_sparse$load_npz(io$distMat)

rna <- readRDS(io$rna)
raw_spat <- readRDS(io$spatLocs)

distmat <- csr_matrix
rownames(distmat) <- raw_spat$cell_ID
colnames(distmat) <- raw_spat$cell_ID

Object.sparse <- GetAssayData(rna, slot = "counts",assay = "RNA")

test <- rbind(Object.sparse, distmat)

temp_SO <- CreateSeuratObject(counts=test)

assayName <- "RNA"
normalizationMethod <- "CLR"
ntopics <- as.numeric(opts$ntopics)
alpha <- 50
beta <- 0.1
varFeatures <- 5000
iterations <- 500
burnin <- 250
seed.number <- 8
outDir <- io$outdir
cores <- 1

Object        <- NormalizeData(temp_SO, assay = "RNA", normalization.method = "CLR")
#Object        <- FindVariableFeatures(Object, assay = assayName, nfeatures = 5000)
Object.sparse <- GetAssayData(Object, slot = "data",assay = "RNA")
#Object.sparse <- Object.sparse[VariableFeatures(Object, assay = assayName),]

  #convert data into the proper input format for lda.collapsed.gibbs.sampler
data.use      <- Matrix::Matrix(Object.sparse, sparse = T)

data.use      <- data.use * 10
data.use      <- round(data.use)
data.use      <- Matrix::Matrix(data.use, sparse = TRUE)
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
     if (!dir.exists(outDir)) {
       dir.create(outDir)
     }
     saveRDS(selected.Model, paste0(outDir, "/Model_", as.character(topics), "topics.rds"))
}

model_maker(ntopics)
