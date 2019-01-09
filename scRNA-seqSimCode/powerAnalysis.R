#Douglas Abrams
#The Jackson Lab, 8.4.17

#executes the expiremental power analysis
#clusterAg = cluster algorithm ("k-means", "sc3", "sincera", "seurat", "snn", "SIMLR", "pcaReduce" )
#normAlg = normalization algorithm ("scran", "scnorm", "seurat", "TMM", "DESeq", "RPM")
#threshold for acceptable true positive rate. suggested:0.8
#cellRarity: rarity of minority cell type
#lower.bound = minimum cell type. Suggested: 1
#upper.bound = maximum cell type. Suggested: 200/ cellRarity ... at max, 200 rare cells
#midpoints: a list to hold the midpoints to detect if the algorithm has converged on a set of sample sizes

source("pipeline.R")
powerAnalysis = function(ClusterAlg, normAlg, FC, nGenes, threshold , cellRarity,  lower.bound, upper.bound =  200/cellRarity, midpoints=c() ){
  
  if (lower.bound < upper.bound){
    #calculate midpoint sample size
    mid = lower.bound + (upper.bound - lower.bound)/2
    
    #if the midpoint is staying the same then return it
    midpoints = append(midpoints, mid)
    if (length(midpoints) > 5 && floor(midpoints[length(midpoints)]) ==floor(midpoints[length(midpoints)-1]) && floor(midpoints[length(midpoints)]) == floor(midpoints[length(midpoints)-1])){
      return(mid)
    }
    #add mid to list
    
    #calculate groupCells param
    groupCells = c( round(mid - (mid * cellRarity),0), round(mid * cellRarity,0) )
    
    #list of Rands
    Rands=unlist(list())
    
    #run pipeline three times
    for (i in 1:3){
      RI = pipeline(ClusterAlg, normAlg ,groupCells, FC, nGenes)
      print(RI)
      Rands[i] = RI$Randindx
    }
    
    #take average Rand
    Rand = mean(Rands)
    
    print(paste0(Rand, " rand with ",mid))

    #recurse etc.
    if (Rand < threshold){
      if (threshold - Rand < 0.01) {
        return(list("Rand" = Rand, "pCells " = mid ))
      }
      else{
        powerAnalysis(ClusterAlg, normAlg, FC, nGenes, threshold, cellRarity, mid, upper.bound, midpoints)
      }
    }
    else{
      if (Rand - threshold < 0.01){
        return(list("Rand" = Rand, "pCells " = mid ))
      }
      else{
        powerAnalysis(ClusterAlg, normAlg, FC, nGenes, threshold, cellRarity,  lower.bound, mid, midpoints)
      }
    }
  }
}






