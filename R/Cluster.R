#' Normalization
#'
#' This function serves as the data preprocessing step for SHARP
#'
#' @param X a count expression matrix.
#'
#' @return A normalized single-cell expression matrix.
#'
#' @examples
#'
#' XNorm <- NormalizeSC(X)
#'
#'
NormalizeSC <- function(X) {
  k <- colSums(X)
  k <- k/median(k)
  Xnorm <- sweep(X, 2, k, FUN = "/")
  return(list(normdata = Xnorm))
}



#' Compute initial value
#'
#' @param Y a count expression matrix without missing value.
#' @param celltype a cell label array, which could be based on other clustering methods or prior biological knowledge.
#' @param subjectid a subject ID array.
#'
#' @return A list containing the initial value of cell-type effect based on each subject.
#' @export
#'
#' @examples
#'
#' SLabel <- SHARP(data, N.cluster = Ncluster, rN.seed = 1234)
#'
#' alpha_0 <- InitVal(data, SLabel$pred_clusters, subjectid)$alpha[[baseline_ID]]
#'
#'
InitVal <- function(Y, celltype, subjectid) {
  ## normalize the counts by total
  k <- colSums(Y)
  k <- k/median(k)
  Ynorm <- sweep(Y, 2, k, FUN = "/")

  ## separate by subject
  subjects <- unique(subjectid)
  nsub <- length(subjects)
  alpha <- vector("list", nsub)

  celltypes <- unique(celltype)
  ncell <- length(celltypes)
  ngene <- nrow(Y)

  for (i in 1:nsub) {
    ix <- subjectid == subjects[i]
    thisY <- Ynorm[, ix]
    thisCell <- celltype[ix]
    alpha[[i]] <- matrix(0, ngene, ncell)
    for (j in 1:ncell) {
      ix.cell <- thisCell == celltypes[j]
      if (sum(ix.cell) > 0) {
        tmp <- rowMeans(thisY[, ix.cell, drop = FALSE])
        alpha[[i]][, j] <- tmp/sum(tmp) * 10
      } else {
        tmp <- rowMeans(thisY)
        alpha[[i]][, j] <- tmp/sum(tmp) * 10
      }
    }
  }
  return(list(alpha = alpha))

}




#' Run SHARP for initial value
#'
#' @param count_all_notna a count expression matrix.
#' @param subject_all_notna a subject ID array.
#' @param Ncluster number of clusters for the final clustering results.
#' @param ID baseline subject ID number.
#' @param seed a number used for setting seeds for SHARP to obtain reproducible results.
#'
#' @return An array containing the initial value of cell-type effect based on baseline subject.
#' @export
#' @import SHARP
#'
#' @examples
#'
#' alpha_0 <- InitVal_S(count_all_notna, subject_all_notna)
#'
#' alpha_0 <- InitVal_S(count_all_notna, subject_all_notna, Ncluster = 5, ID = 3)
#'
#'
InitVal_S <- function(count_all_notna, subject_all_notna, Ncluster = NULL, ID = 1, seed = 1234) {
  count <- count_all_notna[, which(subject_all_notna == levels(factor(subject_all_notna))[ID])]
  subjectid <- subject_all_notna[which(subject_all_notna == levels(factor(subject_all_notna))[ID])]

  data <- NormalizeSC(count)
  data <- log2(data$normdata + 1)

  SLabel <- SHARP(data, N.cluster = Ncluster, rN.seed=seed)

  invisible(capture.output(alpha_0 <- InitVal(count, SLabel$pred_clusters, subjectid)$alpha[[1]]))
  return(alpha_0)
}


#' Run EDClust for single-cell RNA data clustering
#'
#' EDClust: An EM-MM hybrid method for cell clustering in multiple-subject single-cell RNA sequencing
#'
#' @param count_all_notna a count expression matrix without missing value.
#' @param subject_all_notna a subject ID array without missing value.
#' @param alpha_0 the initial value of cell type effect based on InitVal_S() or InitVal()
#' @param EM_steps the maximum number of EM iterations. By default, EM_steps = 100, i.e.,
#' EDClust() will stop and return clustering results when EM_steps = 100.
#' @param MM_steps the maximum number of MM iterations within the M-step in each EM iteration. By default, MM_steps = 3.
#' @param BaseID is same as the ID used in generating alpha_0. If there's only one sample, you should use the default setting: 0L.
#' @param flag the stopping_criterion of EDClust. By default, flag = 3.
#'
#' @return EDClust returns a list object containing:
#' @return \itemize{
#'   \item mem: an array of clustering label
#'   \item loglik: the final log likelihood after iterations
#'   \item alpha_0: a matrix of cell-type specific effects estimates
#'   \item delta: a three dimensions matrix of subject-type specific effects estimates
#'   \item alpha: a three dimensions matrix of overall effects estimates
#'   \item p: a matrix with probability that each cell belongs to each cluster
#' }
#'
#' @export
#'
#' @examples
#'
#' EMMM_Result <- FitPolya(count_all_notna, subject_all_notna, alpha_0)
#'
#' EMMM_Result <- FitPolya(count_all_notna, subject_all_notna, alpha_0,
#'                         EM_steps = 20L, MM_steps = 5L, BaseID = 2L, flag = 1e-05)
#'
#'
FitPolya <- function(count_all_notna, subject_all_notna, alpha_0, EM_steps = 100L, MM_steps = 3L, BaseID = 0L, flag = 1e-04) {
  EMMM_Result <- julia$call("fitPolya", count_all_notna, subject_all_notna, alpha_0, EM_steps, MM_steps, BaseID, flag)
  names(EMMM_Result) <- c("mem", "loglik", "alpha_0", "delta", "alpha", "p")
  return(EMMM_Result)
}







