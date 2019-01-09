library

#generates counts data. k = number of groups. number = number of cels per group. Type = the type of counts data (whether or not
#to include a numeric offset between group populations. Difference = numeric difference
#between number of cells in each group (i.e. if type = "ranged,"number = 100 & difference = 5, group 1 = 95, group 2 = 105). 

generateCounts = function(cellGroups = c(100,100),batch = F){

  if (is.null(cellGroups) == T){
    cellGroups = c(100,100)
  }
  cellGroups = unlist(cellGroups)

  #simulate base data
  if(batch == F){
    splat = newSplatParams(nGenes = 10000, batchCells = sum(cellGroups), dropout.present = T) 
    sim = splatSimulate(splat, verbose = F)
  }
  if (batch ==T){
    splat = newSplatParams(nGenes = 10000, batchCells = cellGroups, dropout.present = T, batch.facLoc =c , batch.facScale=d) 
    sim = splatSimulate(splat, verbose = F)
  }


  real_counts = unlist(list())
  for (i in 1:length(cellGroups)){
    real_counts=append(real_counts, rep(i, cellGroups[i]))
  }
  return(list("sim" = sim, "groupCells" = cellGroups, "splat" =  splat, "clusterGroups" = real_counts))
}


