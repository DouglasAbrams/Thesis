

source("getVars.R")

#clusters data using Seurat. All filtering parameters have been set such that filtering wont occur
seurat = function(counts) {
  if (dim(counts)[1] < dim(counts)[2]){
    counts = t(counts)
  }
  counts = log(counts+1)
  seur <-CreateSeuratObject(raw.data = counts, min.cells = 3, min.genes = 200, project = "project")
  #highly dispersed genes
  seur <- FindVariableGenes(object = seur, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  print(seur@var.genes)
  
  seur <- ScaleData(object = seur)
  

  #PCA & tsne
  seur <- RunPCA(object = seur, pc.genes = seur@var.genes, do.print = T, pcs.print = 1:5, 
                 genes.print = 5)
  #clustering
  
  seur <- JackStraw(object = seur, num.replicate = 100, do.print = FALSE)
  
  seur <- FindClusters(object = seur, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 0.6, print.output = 0, save.SNN = TRUE)
  print(seur@ident)
  return(seur@ident)
}

#Find best parameter G. 20 iterations.

###_____THIS function is no loner used as of Seurat 2.0#############
optomizedRunSEURAT = function(data){
  print("running seurat 20 times at different densities")
  ps = unlist(list())
  Gs = c()
  
  for (G in Gs){
    clust = seurat(data, G)
    print(clust)
    p = getVariance(data, clust)
    ps[G] =  p
  }
  return(seurat(data, G = which.min(ps)))
}
