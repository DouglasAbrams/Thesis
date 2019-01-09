source("C:/Users/s-abramd/Desktop/SSP_Proj_2017/code_for_pipeline/pipeline.R")


#run on either "manyLowDE obj" or "fewhighDE" obj to test clustering methods
#they were generated using the "pipeline" function

#create two datasets - one with many lowly expressed DE genes and one with 
#few lowly expressed ones




for (i in 6:10){
  resultsKM[i] = cluster(fewHighDE,"k-means")
}

for (i in 6:10){
  resultsSN[i] = cluster(fewHighDE,"snn")
}

for (i in 6:10){
  resultsSC[i] = cluster(fewHighDE,"sc3")
}

for (i in 6:10){
  resultspcaR[i] = cluster(fewHighDE,"pcaReduce")
}

for (i in 6:10){
  resultsSIM[i] = cluster(fewHighDE,"SIMLR")
}

for (i in 6:10){
  resultsSINC[i] = cluster(fewHighDE,"sincera")
}

for (i in 6:10){
  resultsSeur[i] = cluster(fewHighDE,"seurat")
}
