
#runs k-means
#reduced_counts  = counts matrix post PCA
#k = numClusters
kmeans_cluster = function(reduced_counts, k){
  means_clusters = stats::kmeans(reduced_counts, k, nstart = 1)
  centers = means_clusters$centers
  clust=means_clusters$cluster 

  return(clust)
}

?kmeans
