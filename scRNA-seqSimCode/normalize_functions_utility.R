# gene filtering

"****directly taken to test from normalizing single cell RNA.... 
Catalina Vallejos et. al." 

#Functions #####
#RPM scaling factor
rpm <- function(x) {
  print(dim(x))
  s <-colSums(x)
  n.cells<-ncol(x)
  s <- n.cells*(s/sum(s))
  counts <-t(t(x)/s)
  print(dim(counts))
  return(list(s=s, counts=counts))
}

#deseq normalization
deseq <- function(x) {
  print(dim(x))
  
  s <-estimateSizeFactorsForMatrix(x)
  n.cells<-ncol(x)
  s <- n.cells*(s/sum(s))
  counts <-t(t(x)/s)
  print(dim(counts))
  
  return(list(s=s, counts=counts))
}

tmm <- function(x) {
  s <-edgeR::calcNormFactors(x)*colSums(x)
  n.cells<-ncol(x)
  s <- n.cells*(s/sum(s))
  counts <-t(t(x)/s)
  return(list(s=s, counts=counts))
}

