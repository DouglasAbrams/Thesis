#code for performing the various tests on normalization methods


source("generateCounts.R")
source("normalize.R")
source("DE5.R")
source("QCfilter.R")
source("normalizeTestDataSets.R")
source("addNoise")
source("testNormalizationGroupings")




#function that tests normalizations with MAST
#BASE = obj that contains data + simulation objects
#do compare = T when you want TPR F when you want raw value of num genes detected
testNorm = function(BASE, normalize = T,  normAlg, doCompare = T){              
  
  if (normalize == T){
    print("***NORMALIZING***")
    counts = normalize(normAlg, BASE) 
  }
  return(assessNorm(counts, BASE$sim, BASE$markerGenes, doCompare))
}

#takers in a BASE obj with inflated library sizes and applied a normalization method,
#calculating the percent to which the method lowered the lib. size variation
#normAlg = algorithm for normalization
#BASE = base object
#nonNormedObj = object before normalization
reduceSkew = function(normAlg, BASE, nonNormedOBJ){
  counts = normalize(normAlg, BASE) 
  
  sums = colSums(counts)
  sumsList = c()
  for (k in 1:length(sums)){
    sumsList[k] = sums[k]
  }
  print(length(sumsList))
  print(dim(counts))
  normalizedGroupings = as.data.frame(matrix(nrow = dim(counts)[2], ncol = 2))
  normalizedGroupings[,1] = sumsList
  print("here")
  
  normalizedGroupings[,2] = nonNormedOBJ$cellValsFrame$group[!BASE$outlierCells]
  colnames(normalizedGroupings) = c("colSums", "group")
  
  varianceChange = testNormalizationGroups(nonNormedOBJ, normalizedGroupings)
  return(varianceChange)
}

#parameters for generating test datasets
nCells_pGroup = c(200,200)
nGenes_pGroup = c(50,50)
FC_pGroup = c(2,-2)
levels = c(1, 4)


#code for generating sample datasets for tests
print("***GENERATING INITIAL SCEset WITH COUNTS***")
initial = generateCounts(nCells_pGroup,F)
print("***GENERATING SPECIFIC DE GENES***")

#use to test marker genes
BASE = createdeGENES_V.3(initial$sim, initial$splat, nGenes_pGroup, FC_pGroup, nCells_pGroup)

#use for first test
output = addNoise(BASE, levels)
BASE = output$BASE

#for exact change
unlockBinding( "counts", assayData(initial$sim))
assayData(initial$sim)$counts[, 201:400] = assayData(initial$sim)$counts[,201:400] *4
BASE = initial

pData(BASE$sim)$condition = rep(1:length(nCells_pGroup),nCells_pGroup)

markerGeneInformation = BASE$information_matrix

print("***FILTERING AND RUNNING QC***")
BASE = QCFilter(BASE$sim,markerGeneInformation)



####test for no marker genes###################################

pData(initial$sim)$condition = rep(1:length(nCells_pGroup),nCells_pGroup)

markerGeneInformation = initial$information_matrix

print("***FILTERING AND RUNNING QC***")
BASE = QCFilter(initial$sim,markerGeneInformation)

#test to see how many marker genes MAST detects
noMG = c()
noMG[1] = nrow(testNorm(BASE, T, "scran", doCompare = F))
noMG[2] = nrow(testNorm(BASE, T, "DESeq", doCompare = F))
noMG[3] = nrow(testNorm(BASE, T, "RPM", doCompare = F))
noMG[4] = nrow(testNorm(BASE, T, "TMM", doCompare = F))
noMG[5] = nrow(testNorm(BASE, T, "Seurat", doCompare = F))
noMG[6] = nrow(testNorm(BASE, T, "scnorm", doCompare = F))



##########################tests for 50 marker genes + detect 1/ MAST##############

#here we test four datasets, with 0.5 MG FC, 1, 1.5 and 2

#test Seurat
seur = c()
seur[1] = testNorm(p5, T, "seurat")$TPR
seur[2] = testNorm(one, T, "seurat")$TPR
seur[3] = testNorm(onep5, T, "seurat")$TPR
seur[4] = testNorm(two, T, "seurat")$TPR

TMM = c()
TMM[1] = testNorm(p5, T, "TMM")$TPR
TMM[2] = testNorm(one, T, "TMM")$TPR
TMM[3] = testNorm(onep5, T, "TMM")$TPR
TMM[4] = testNorm(two, T, "TMM")$TPR

RPM = c()
RPM[1] = testNorm(p5, T, "RPM")$TPR
RPM[2] = testNorm(one, T, "RPM")$TPR
RPM[3] = testNorm(onep5, T, "RPM")$TPR
RPM[4] = testNorm(two, T, "RPM")$TPR

DESeq = c()
DESeq[1] = testNorm(p5, T, "DESeq")$TPR
DESeq[2] = testNorm(one, T, "DESeq")$TPR
DESeq[3] = testNorm(onep5, T, "DESeq")$TPR
DESeq[4] = testNorm(two, T, "DESeq")$TPR

scran = c()
scran[1] = testNorm(p5, T, "scran")$TPR
scran[2] = testNorm(one, T, "scran")$TPR
scran[3] = testNorm(onep5, T, "scran")$TPR
scran[4] = testNorm(two, T, "scran")$TPR

scnorm = c()
scnorm[1] = testNorm(p5, T, "scnorm")$TPR
scnorm[2] = testNorm(one, T, "scnorm")$TPR
scnorm[3] = testNorm(onep5, T, "scnorm")$TPR
scnorm[4] = testNorm(two, T, "scnorm")$TPR


##########################tests for % reduced variation##############



#run reduceSkew

change1 = reduceSkew("scran", BASE, output)
change2 = reduceSkew("DESeq", BASE, output)
change3 = reduceSkew("RPM", BASE, output)
change4 = reduceSkew("TMM", BASE, output)
change5 = reduceSkew("Seurat", BASE, output)
change6 = reduceSkew("scnorm", BASE, output)

changeVARS[7] = change1
changeVARS[8] = change2
changeVARS[9] = change3
changeVARS[10] = change4
changeVARS[11] = change5
changeVARS[12] = change6

comparison = as.data.frame(matrix(nrow = 6, ncol = 2))
comparison[,1] = changeVARS[1:6]
comparison[,2] = changeVARS[7:12]
rownames(comparison)  = c("scran", "DESeq", "RPM", "TMM", "Seurat", "SCnorm")
colnames(comparison) = c("4 groups", "8 groups")