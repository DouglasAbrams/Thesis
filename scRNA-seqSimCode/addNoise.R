#inflates the library sizes of a subset of cell groups (# of cell groups and magnitude of scaling is specified by "levels")
#BASE = base object with data sets and sim objects
#levels = specifid like other _pGroup functions. So if I want 3 groups with 1,4, and 8 lib.
#size inflation I would input c(1,4,8)
addNoise = function(BASE,levels){
  
  #get variables ready
  sim = BASE$sim
  counts = assayData(sim)$counts
  orC = counts
  nGroups = length(levels)
  
  #assemble list of group/group inflations
  groupSize = floor(dim(counts)[2]/nGroups)
  remainder = dim(counts)[2] - (groupSize * nGroups)
  
  breaks = rep(groupSize, nGroups)
  breaks[1] = breaks[1] + remainder
  skewLEVS = rep(0,dim(counts)[2])
  
  #skew library sizes w/ negative binomial
  for (i in 1:length(breaks)){

    min = sum(breaks[1:i-1]) + 1
    max = breaks[i] + sum(breaks[1:i-1])
    print(min)
    print(max)
    lev = levels[i]
    if (lev < 0){
     lev = lev *-1
     skews = rnbinom(breaks[i], mu = lev, size = 1)
     skews = skews *-1
     for (j in 1:length(skews)){
       if (skews[j] == 0){
         skews[j] = -1
       }
     }
    }
    else{
      skews = rnbinom(breaks[i], mu = lev, size = 1)
    }
    for (j in 1:length(skews)){
      if (skews[j] == 0){
        skews[j] = 1
      }
    }
    skewLEVS[min:max] = skews
    counts[,min:max] = t(sweep(t(counts[,min:max]), MARGIN = T, skews, '*'))
  }
  sums = colSums(counts)
  sumsList = c()
  for (k in 1:length(sums)){
    sumsList[k] = sums[k]
  }
  unlockBinding("counts", assayData(BASE$sim))
  assayData(BASE$sim)$counts = counts
  
  dataF = as.data.frame(matrix(nrow = sum(breaks), ncol = 2))
  dataF[,1] = sumsList
  dataF[,2] = rep(c(1:length(breaks)),breaks)
  colnames(dataF) = c("colSums", "group")
                    
  return(list("BASE"=BASE,"FCs" = skewLEVS, "cellValsFrame" = dataF))
}

