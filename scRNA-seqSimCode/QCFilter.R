#QC and filters a SCeset object using calculateQCMetrics (splater). It checks for library size & coverage. 

QCFilter = function(sim,markerGenes){
  
  #filter genes
  ave.counts <- calcAverage(sim)
  keep <- ave.counts >= 1
  
  #dont filter marker genes
  for (index in markerGenes$gene_index){
    keep[index] = TRUE
  }
  #filter cells
  #calculate QC on both depth and coverage, outlier >= 3 MADs away
  QC = calculateQCMetrics(sim)
  
  librarySizeOutliers = isOutlier(QC$total_counts, nmads = 3, log = TRUE)
  coverageOutliers = isOutlier(QC$total_features, nmads = 3, log = TRUE)
  
  for (i in 1:length(coverageOutliers)){
    if (coverageOutliers[i] != librarySizeOutliers[i] && coverageOutliers[i] == FALSE){
      coverageOutliers[i] = librarySizeOutliers[i]
    }
  }
  outliers = coverageOutliers
  
  #set simulation
  sim = sim[keep,!outliers]
  
  #update marker gene indexes based on filtering
  for(i in 1:length(fData(sim)[,1])){
    if (as.character(fData(sim)[,1][i]) %in% markerGenes$gene){
      markerGenes$gene_index[markerGenes$gene == as.character(fData(sim)[,1][i])] = i
    }
  }
 return(list("sim" = sim, "markerGenes" = markerGenes, "outlierCells" = outliers))
}
