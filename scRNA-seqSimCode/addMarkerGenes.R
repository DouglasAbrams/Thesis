#this function adds specific marker genes at specific fold changes to a counts matrix
#the value added is the fold change * or / (depending on + or -) to the median of all
#counts values that are not in the group being changed. 

#sim --> a splat simulation
#splat --> splat params

##for each of these a group is specified as an index in an atomic. 

##i.e. FC_pGroup = (4, -6, 0.8, 2) would mean that cell group one has a FC of 4 and so on. 

#num_Mark_Genes_pGroup--> the number of DE genes to create in each cell type population
#FC_pGrooup -->> the center of a distribution from which to sample FC values for each cell type population
#groupCells --> the number of cells in each cell type group




createdeGENES_V.3 = function(sim,splat, num_Mark_Genes_pGroup, FC_pGroup, groupCells){
  
  gene_indexes_pGroup = unlist(list())
  infoMatrix = matrix(nrow = sum(num_Mark_Genes_pGroup), ncol = 6)
  data = assayData(sim)$TrueCounts

  matrix_withMedians = getIndexMedians(num_Mark_Genes_pGroup, infoMatrix, groupCells, data)
  matrix_withFCs = simulateFold_changes(matrix_withMedians, FC_pGroup)
  res = addToMatrix(groupCells, matrix_withFCs, data,sim)
  
  data_with_FC = res$data
  infoMatrix = res$infoMatrix
  infoMatrix = as.data.frame(infoMatrix)
  colnames(infoMatrix) = c("gene_index","median", "group", "fold_change", "trueCount_value","gene")
  infoMatrix$gene = paste0("Gene", as.character(infoMatrix$gene_index))
  print(infoMatrix)
  set_exprs(sim, "TrueCounts") <- data_with_FC
  final_sim = splatSimDropout(sim, splat)
  
  return(list("sim" = final_sim, "information_matrix" = infoMatrix))
  
}

getIndexMedians = function(num_Mark_Genes_pGroup, infoMatrix, groupCells, original_means) {
  if (length(num_Mark_Genes_pGroup) == 1) {
    print("need at least two groups")
    return
  }
  rm(.Random.seed, envir=globalenv())
  
  nGenes = sum(num_Mark_Genes_pGroup)
  indexes = unlist(round(runif(nGenes, min = 1, max = dim(original_means)[1]), 0 ))
  
  return(assembleMatrix(indexes, num_Mark_Genes_pGroup, groupCells, original_means, infoMatrix))
}

assembleMatrix =  function(indexes, num_Mark_Genes_pGroup, groupCells, original_means, infoMatrix){
  tempMarker = append(0,num_Mark_Genes_pGroup)
  tempGroup = append(0, groupCells)
  for (k in 2:length(tempMarker)){
    min = sum(tempMarker[1:k-1]) + 1
    max = tempMarker[k] + sum(tempMarker[1:k-1])
    
    genes_indexes = indexes[min:max]
    for ( t in 1:length(genes_indexes)){a
      index =genes_indexes[t]
      min = sum(tempGroup[1:k-1]) + 1
      max = tempGroup[k] + sum(tempGroup[1:k-1])
      temp = unlist(list())
      for (i in min:max) {
        temp[i] = i
      }
      median = mean(original_means[index, -temp[min:max]])
      if(median==0){
        median = 1
      }
      infoMatrix[t+sum(tempMarker[1:k-1]), 1] = index    
      infoMatrix[t+sum(tempMarker[1:k-1]), 2] = median    
      infoMatrix[t+sum(tempMarker[1:k-1]), 3] = k - 1  
    }
  }
  return(infoMatrix)
}

simulateFold_changes <-function(infoMatrix,FC_pGroup) {
  
  for (index in 1:length(FC_pGroup)){
    
    fc = FC_pGroup[index]
    matrixSectionToChange = infoMatrix[infoMatrix[, 3] == index,]
    if (fc < 0){
      fc = fc * -1
      foldChanges = rnbinom(length(matrixSectionToChange[,4]), mu = fc, size = 1)
      foldChanges = -1 * foldChanges
      fc = fc * -1 
    }
    if (fc >= 0){
      foldChanges = rnbinom(length(matrixSectionToChange[,4]), mu = fc, size = 1)
    }

    infoMatrix[infoMatrix[, 3] == index,4] = foldChanges

  }
  return(infoMatrix)
}

#appends to the matrix the DE genes i..j Einf --> log(i/Jsoid)
addToMatrix <- function(groupCells,infoMatrix,original_means,sim) {
  
  tempList = list()
  for (m in 1:length(groupCells)){
    
    min = sum(groupCells[1:m-1]) + 1
    max = groupCells[m] + sum(groupCells[1:m-1])
    
    indexes  = infoMatrix[infoMatrix[,3]== m ,1]
    fold_changes  = infoMatrix[infoMatrix[,3]== m ,4]
    medians  = infoMatrix[infoMatrix[,3]== m ,2]
  
    for (i in 1:length(indexes)) {
      
      index = indexes[i]
      FC = fold_changes[i]
      median = medians[i]
      
      if (FC >= 0) {
        original_means[index, min:max]  =  median * FC
        original_means[index, min:max] = original_means[index, min:max] *(pData(sim)$ExpLibSize[min:max]/mean(pData(sim)$ExpLibSize[min:max]))
        tempList = append(tempList,median * FC)
      }
      if (FC < 0) {
        original_means[index, min:max]  =  median / (-1 * FC)
        original_means[index, min:max] = original_means[index, min:max] *(pData(sim)$ExpLibSize[min:max]/mean(pData(sim)$ExpLibSize[min:max]))
        tempList = append(tempList,median / (-1 * FC))
        
      }
      
    }
  }
  infoMatrix[, 5] =  unlist(tempList)
  return(list("data" = original_means, "infoMatrix" = infoMatrix))
}

splatSimDropout <- function(sim, params) {
  dropout.present = F
  true.counts <- get_exprs(sim, "TrueCounts")
  
  if (dropout.present) {
    cell.names <- pData(sim)$Cell
    gene.names <- fData(sim)$Gene
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    dropout.mid <- getParam(params, "dropout.mid")
    dropout.shape <- getParam(params, "dropout.shape")
    cell.means <- get_exprs(sim, "CellMeans")
    
    # Generate probabilites based on expression
    drop.prob <- sapply(seq_len(nCells), function(idx) {
      eta <- log(cell.means[, idx])
      return(logistic(eta, x0 = dropout.mid, k = dropout.shape))
    })
    
    # Decide which counts to keep
    keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob),
                   nrow = nGenes,
                   ncol = nCells)
    
    counts <- true.counts * keep
    
    colnames(drop.prob) <- cell.names
    rownames(drop.prob) <- gene.names
    colnames(keep) <- cell.names
    rownames(keep) <- gene.names
    
    set_exprs(sim, "DropProb") <- drop.prob
    set_exprs(sim, "Dropout") <- !keep
  } else {
    counts <- true.counts
  }
  
  scater::counts(sim) <- counts
  
  return(sim)
}

logistic <- function(x, x0, k) {
  1 / (1 + exp(-k * (x - x0)))
}
