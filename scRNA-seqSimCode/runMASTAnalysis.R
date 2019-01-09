#finds DE genes in expression matrix counts from object sim
#markerGenes = marker gene matrix
#doCompare = whether or not to return raw num. DE genes or TPR

assessNorm = function(counts, sim, markerGenes, doCompare){
  print(dim(counts))
  

  #find DE genes with MAST
  
  #construct object
  mastOBJ = FromMatrix(counts,pData(sim), fData(sim) )
  
  #find the # genes detected in a sample
  cond = factor(colData(mastOBJ)$condition)
  cdr2 = colSums(assay(mastOBJ)>0)
  colData(mastOBJ)$cngeneson <- cdr2
  
  #run a zero-inflated regression
  zlmCond <- zlm(~condition + cngeneson, mastOBJ)
  
  #compile DE genes, ranked by descrete z score, correlation and a continous z score based on model
  #pull a list the length of the number of marker genes
  summaryCond <- summary(zlmCond, doLRT='condition') 
  
  
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[contrast=='condition' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='condition' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle[fdr<.05], as.data.frame(mcols(mastOBJ)), by='primerid')
  
  colnames(fcHurdleSig)[1] = "gene"

  if (doCompare == T){
  TPR = length((fcHurdleSig$gene[fcHurdleSig$gene %in% markerGenes$gene == T]))/ length(fcHurdleSig$gene) 
  return(list("TPR" = TPR, "FPR" = 1-TPR))
  }
  return("MGMat"= fcHurdleSig)
}
