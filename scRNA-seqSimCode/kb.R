#performs a bubble sort on any data type with order specified by user
#order = 1, ascending
#order = 0, descending
bubbleSort = function(elems, order){
  if (order != "ascending" && order != "descending"){
    print("Order Usage:<order> = 'ascending' or 'descending'")
    
  }
  sortCount = 1
  if (order == 1){
    while (sortCount > 0){
      sortCount = 0
      for (i in 1:(length(elems) - 1)){
        if (elems[i] > elems[i+1]){
          first = elems[i]
          second = elems[i+1]
          elems[i] = second
          elems[i+1] = first
          sortCount = sortCount + 1;
        }
      }
    }
    return(elems)
  }
  if (order == 0){
    while (sortCount > 0){
      sortCount = 0
      for (i in 1:(length(elems) - 1)){
        if (elems[i] < elems[i+1]){
          first = elems[i]
          second = elems[i+1]
          elems[i] = second
          elems[i+1] = first
          sortCount = sortCount + 1;
        }
      }
    }
    return(elems)
  }
}

t = c(c(1,2,2,2), c(10000000,3,5,6), c(375,667))
new = bubbleSort(t,1)


#
