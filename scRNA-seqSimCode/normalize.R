
#normalizing QC obj, which contains all preceeding information matrices and sim objects
#alg = normalization algorithm
source("normalize_functions_utility.R")
normalize = function(alg,QCobj){
  sim = QCobj$sim
  counts = counts(sim)
  
  #unlock sceset binding
  unlockBinding("counts", assayData(sim)) 
  
  if(alg == "scran"){
    if (dim(counts(sim))[2] > 500 || length(nCells_pGroup) >= 4){
      clusters <- quickCluster(counts)
      counts <- computeSumFactors(counts, cluster=clusters)
    }
    else{
      factors = computeSumFactors(sim) 
      sce <- scater::normalise(factors)
      counts = assayData(sce)$norm_exprs
     # counts = t(t(counts)/factors*mean(factors)) 
    }

  }
  if(alg == "scnorm"){
    conditions = pData(sim)$condition
    counts = SCnorm(Data = counts, Conditions = conditions, OutputName = "normedData", FilterCellNum = 10, PropToUse = 0.25, Thresh = 1, ditherCounts = F, K = length(nCells_pGroup))$NormalizedData

  }
  if(alg == "RPM"){
    if (dim(counts)[1] <  dim(counts)[2]){
      counts = t(counts)
    }
    counts = rpm(counts)$counts
  }
  if(alg == "DESeq"){
    if (dim(counts)[1] <  dim(counts)[2]){
      counts = t(counts)
    }
    counts = deseq(counts)$counts

  }
  if(alg == "TMM"){
    if (dim(counts)[1] <  dim(counts)[2]){
      counts = t(counts)
    }
    res = tmm(counts)
    counts = res$counts
    print(res$s)
  }
  
  if (alg == "seurat"){
    #make data into seurat obj
    seur <- CreateSeuratObject(raw.data = counts, project = "sampleSeurat")
    seur <- NormalizeData(object = seur, normalization.method = "LogNormalize", scale.factor = 10000)
    counts = as.matrix(seur@data)
  }
  print(dim(counts))
  return(counts)
}
