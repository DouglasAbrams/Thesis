
#executes a sincera cluster. the parts that are important are the parameters used for selecting genes. (specificity.thresh)

sincera_function = function(counts, real_counts, k){

  #data orientation
  if (dim(counts) [1] <  dim(counts)[2]){
    counts = t(counts)
  }    
  print(dim(counts))
  print(length(real_counts))
  counts = counts*counts
  sincera_obj <- construct(exprmatrix =as.data.frame(counts), samplevector=paste("sample", real_counts, sep=""))
  sincera_obj <- normalization.zscore(sincera_obj, pergroup=TRUE)
  
  #select genes
  sincera_obj <- cluster.geneSelection(sincera_obj, method="specificity", pergroup=TRUE, min.samples=2, specifity.thresh=0.7)
  
  #dim reduction
  sc <- doPCA(sincera_obj, genes=getGenesForClustering(sincera_obj), use.fast = T, center=T, scale.=T)
 # sc <- doTSNE(sc, dims=1:5)
  
  sincera_obj = cluster.assignment(sc,k=k)
  clust = getCellMeta(sincera_obj, "GROUP")  
  
  return(clust)
}
