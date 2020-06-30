#' GeneScoreDistubiton
#'
#'
#' This functionreturns a matrix of genescores for each topic
#'
#' @param Model LDA Model
#'
#' @examples
#' GeneScores(Model)
#'
#' @return Vmatrix of Gene-Topic Scores
#'
#'
#' @export
#'
#' @import RColorBrewer


GeneScores <- function(Model) {
  normalized.topics <- Model$topics/(rowSums(Model$topics) + 1e-05)
  scores <- apply(normalized.topics, 2, function(x) x *
                    (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
  scores <- t(scores)
  colnames(scores) <- paste("Topic", 1:ncol(scores), sep = "_")
  return(scores)
}
