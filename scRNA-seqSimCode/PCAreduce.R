#executes pcaReduce
#This file, unlike the other clustering methods, returns a RI and a confusion matrix,
# in addition to a set of clusterings

#counts = counts matrix

#k = num clusters
#real_Clust = set of true clusterings (used to calculate RI)

pcaReduce_function = function(counts,k, real_clust){
  reduced = PCAreduce(counts,100, k*2, "M")
  width = dim(reduced[[1]])[2]
  maxRam = -1
  maxConf = 0
  clust = c()
  bestClust = c()
  print(dim(counts))
  reduced_counts = prcomp(x = counts)$x
  for (i in 1:length(reduced)){
    output = unlist(reduced[[i]][,width])
    
    for (j in 1:length(output)){
      clust[j] = output[j]
    }
    RI = rand.index(clust, unlist(real_clust))
    confMat = confusionMatrix(clust, unlist(real_clust))
    if (RI>maxRam){
      maxRam = RI
      maxConf = confMat
      bestClust = clust
    }
  }
  plot(reduced_counts, col = bestClust)
  return(list("Randindx" = maxRam, "confusionMatrix" = maxConf, "clust" = bestClust))  
}
