
#execute optomizedRunSNN to run SNN-cliq. Takes in normalized data matrix.
source("getVars.R")

SNN<-function(data, outfile, k, distance){
  
  if(missing(data)){
    stop(paste("Input data missing.",help,sep="\n"))
  }
  if(missing(outfile)){
    stop(paste("Output file name missing.",help,sep="\n"))
  }
  if(missing(k)){
    k=3
  }
  if(missing(distance)){
    distance<-"euclidean"  # other distance options refer to dist() in R
  }
  m<-as.data.frame(data)
  numSpl<-dim(data)[1]
  m<-dist(data, distance)
  x<-as.matrix(m)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  
  edges<-list()              # SNN graph
  for (i in 1:numSpl){
    j<-i
    while (j<numSpl){
      j<-j+1
      shared<-intersect(IDX[i,], IDX[j,])
      if(length(shared)>0){			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
        strength<-max(s)
        if (strength>0)
          edges<-rbind(edges, c(i,j,strength))
      }				
    }
  }
  write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}

snnCliq<-function(par.k=30, par.r=0.8, par.m=0.5, dat, infile, outfile, distance = "euclidean"){
  print(par.m)
  print(par.r)
  d = dat
  SNN(d, "snn-clique.txt", k=par.k, distance)
  system(paste0("python ","C:/Users/s-abramd/Desktop/final_files/clusteringalgorithms/snn_Stuff/",  "Cliq.py -i ", "snn-clique.txt", "  -o ", "res-snn-clique.txt -r ",
                par.r, " -m ", par.m))
  clusts <- read.table(paste0("res-snn-clique.txt"))
  clust = unlist(list())
  for (i in 1:length(clusts)){
    clust[i] = clusts[i]
  }
  return(clusts$V1)
}

optomizedRunSNN =function(data){
  ps = unlist(list())
  par.ms = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  par.rs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  for (r in par.rs){
    for (m in par.ms){
        print(paste0("iteration", m *10))
        clust = snnCliq(par.r = r, par.m = m, dat = data, infile = "snn-cliq.txt", outfile = "res-snn-cliq.txt" )
        print(clust)
        
        p = getVariance(data, clust)
        ps[(m*10) + ((1-(r*10)) * 9)] =  p
      }
  }
  clust = snnCliq(par.m = (which.min(ps)/10), dat = data, infile = "snn-cliq.txt", outfile = "res-snn-cliq.txt" )
  return(clust)
}

