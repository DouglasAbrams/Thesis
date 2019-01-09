#runs the pipeline which filters, normalizes and clusters the data.
#cluster_alg = the clustering algorithm to be used in the analysis
#normAlg = normalization algorithm to be used in the analysis
#nCells_pGroup = the number of cells in each group that is simulated
#FC_pGroup = FC's of the cell groups @ marker genes
#nGenes_pGroup = number of marker genes in each group
#vis = whether or not to plot clustering in PCA space (T/F)

#load files
source("addMarkerGenes.R")
source("QCFilter.R")
source("normalize.R")
source("generateCounts.R")
source("PCAreduce.R")
source("sincera.R")
source("SIMLR.R")
source("k-means.R")
source("sc3.R")
source("Seurat_cluster.R")
source("snn.R")


pipeline = function( cluster_alg, normAlg,  nCells_pGroup =c(600,600), FC_pGroup, nGenes_pGroup){              
  print("***GENERATING INITIAL SCEset WITH COUNTS***")
  initial = generateCounts(nCells_pGroup,F)
  print("***GENERATING SPECIFIC DE GENES***")

  #use to test marker genes
  BASE = createdeGENES_V.3(initial$sim, initial$splat, nGenes_pGroup, FC_pGroup, nCells_pGroup)

  pData(BASE$sim)$condition = rep(1:length(nCells_pGroup),nCells_pGroup)
  
  markerGeneInformation = BASE$information_matrix
  
  print("***FILTERING AND RUNNING QC***")
  QCBASE = QCFilter(BASE$sim,markerGeneInformation)
  
  initial$clusterGroups = initial$clusterGroups[!QCBASE$outlierCells]
  
  unlockBinding("counts", assayData(BASE$sim))
  assayData(BASE$sim)$counts = assayData(BASE$sim)$counts[,!QCBASE$outlierCells]
  
  counts = normalize(normAlg, QCBASE) 
  print("***RUNNING PCA***")
  if (dim(counts)[1] > dim(counts)[2] ){
    counts = t(counts)
  }
  reduced_counts = prcomp(x = counts)$x
  
  infoObj = list()
  infoObj$counts = counts
  infoObj$reduced_counts = reduced_counts
  infoObj$BASE = BASE
  infoObj$initial = initial
  infoObj$nGenes_pGroup = nGenes_pGroup
  infoObj$nCells_pGroup = nCells_pGroup
  
  return(cluster(infoObj, cluster_alg))
}

cluster = function(obj, cluster_alg,vis=T){
  counts = obj$counts
  reduced_counts = obj$reduced_counts
  BASE = obj$BASE
  initial = obj$initial
  nGenes_pGroup = obj$nGenes_pGroup
  nCells_pGroup = obj$nCells_pGroup
  
  if (cluster_alg != "sincera"){
    if (dim(counts) [1] < dim(counts)[2]){
      counts = t(counts)
    }    
  }
  print("***CLUSTERING***")
  if (cluster_alg == "snn"){
    print("***with SNN cliq***")
    clust = optomizedRunSNN(dat = reduced_counts[,1:2] )
    
  }

  if (cluster_alg == "pcaReduce"){
    print ("***with PCA reduce***")
    return(pcaReduce_function(t(counts),length(nCells_pGroup), initial$clusterGroups))
  }

  if (cluster_alg == "sincera"){
    print("***with sincera***")
  #  return(counts)
    clust = sincera_function(counts, initial$clusterGroups, length(nCells_pGroup))
  }
  
  if (cluster_alg == "seurat"){
    print("***WITH SEURAT***")
    clust = seurat(counts)
  }
  
  if (cluster_alg == "SIMLR"){
    print("***WITH SIMLR***")
    counts = t(counts)
    clust = SIMLR(counts,length(nCells_pGroup))
  }
  
  if (cluster_alg == "k-means"){
    print("***WITH K-MEANS***")
    clust = kmeans_cluster(reduced_counts[,1:2],length(nGenes_pGroup))
  }

  if (cluster_alg == "sc3"){
    print("***WITH SC3***")
    ###SC3 takes unfiltered.normalized data because it does it iternally and will error out##
    clust = sc3_function(counts(BASE$sim),length(nCells_pGroup))
  }
  print(dim(reduced_counts))
  if (vis == T){
    plot(reduced_counts, col = clust, xlab = "
         wfnk")
  }
  m = c()
  for (i in 1:length(clust)){
    m[i] = clust[i]
  }
  return(assess_cluster(initial$clusterGroups, m)$Randindex)
}



