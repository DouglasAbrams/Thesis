#Douglas Abrams
#The Jackson Lab, 8.4.17

#assess two cluster assignment. 
#real_counts: simulated assingments
#clust: clustered assignments
#calculates a rand index and confusion matrix

assess_cluster = function (real_counts, clust){
  print(length(real_counts))
  print(length(clust))
  idx = rand.index(real_counts, clust)
  confMat1 = confusionMatrix(clust,real_counts)
  
  for (i in 1:length(clust)){
    if (clust[i] == 1){
      clust[i] = 2
    }
    else{
      clust[i] = 1
    }
  }
  
  confMat2 = confusionMatrix(clust,real_counts)
  
  if (confMat2$overall[1] > confMat1$overall[1]){
    return(list("Randindx" = idx, "confusionMatrix" = confMat2))  
  }
  return(list("Randindx" = idx, "confusionMatrix" = confMat1))  
}
