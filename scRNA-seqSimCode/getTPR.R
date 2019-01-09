#calculates a TPR given an confusion matrix
#confMatriux = confusion matrix
#index = index at which to calculate TPR, set to numRows of matrix

getTPR = function(confMatrix, index){
  print("here")
  return(confMatrix$table[index,index]/ (confMatrix$table[index,index] + confMatrix$table[index, index-1]))
}