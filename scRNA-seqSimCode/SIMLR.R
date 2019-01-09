SIMLR = function(counts, K){
  #load package
  set.seed(11111)
  counts = t(counts)
  clustered = SIMLR_Large_Scale(X = counts,c= K, k=10)
  reduced_counts = clustered$ydata
  clust = clustered$y$cluster
  return(clust)
}


?SIMLR_Large_Scale
