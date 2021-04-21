library(SHARP)
library(mclust)

##Normalization for SHARP
normalizeSC.total <- function(X) {
  k=colSums(X)
  k = k/median(k)
  Xnorm = sweep(X, 2, k, FUN="/")
  list(normdata=Xnorm)
}

#Compute initial alpha0
computeInitVal <- function(Y, celltype, subjectid) {
  ## normalize the counts by total
  k = colSums(Y)
  k = k / median(k)
  Ynorm = sweep(Y, 2, k, FUN="/")
  
  ## seperate by subject
  subjects = unique(subjectid)
  nsub = length(subjects)
  #delta = 
  alpha = vector("list", nsub)
  
  celltypes = unique(celltype)
  ncell = length(celltypes)
  ngene = nrow(Y)
  
  for(i in 1:nsub) {
    ix = subjectid == subjects[i]
    thisY = Ynorm[, ix]
    thisCell = celltype[ix]
    #delta[[i]] = 
    alpha[[i]] = matrix(0, ngene, ncell)
    for(j in 1:ncell) {
      ix.cell = thisCell == celltypes[j]
      if(sum(ix.cell) > 0) {
        tmp = rowMeans(thisY[,ix.cell,drop=FALSE])
        alpha[[i]][,j] = tmp / sum(tmp) * 10
      } else { ## no such cell type in this individual
        ## use the overall mean
        tmp = rowMeans(thisY)
        alpha[[i]][,j] = tmp / sum(tmp) * 10
      }
    }
    #delta[[i]] = matrix(runif(ngene*ncell), ncol=ncell)
    #delta[[i]] = 0.5*alpha[[i]]+1e-10
  }
  
  ## compute initial values for delta. This needs some thinking
  # for(i in 2:nsub) {
  #     tmpDelta = alpha[[i]] - alpha[[i-1]]
  # }
  
  
  #return( list(alpha=alpha, delta=delta) )
  return(list(alpha=alpha))
}
