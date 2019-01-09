#executes an SC3 clustering

#reduced_couynts = counts matrix post PCA
#k = num groups
sc3_function = function(reduced_counts, k){
  #load package

  #orient data
  if(dim(reduced_counts)[1] < dim(reduced_counts)[2]){
    reduced_counts = t(reduced_counts)
  }
  
  #initialize obj
  sce = newSCESet(countData = reduced_counts)

  sc3_obj=calculateQCMetrics(sce)
  
  #run sc3 clustering
  means_cluster = SC3::sc3(sc3_obj, ks = k)
  sc3_clust = pData(means_cluster)[length(pData(means_cluster))]
  clust = list()
  clust = unlist(clust)
  for (i in 1:length(sc3_clust)){
    clust[i] = sc3_clust[i]
  }
  for (a in clust){
    return(a)
  }
}
