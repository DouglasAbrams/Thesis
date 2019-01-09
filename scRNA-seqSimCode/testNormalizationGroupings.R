#takes in two objects, normalized and unnormalizewd and returns 
#the % variance removed via normalization

#unnormalized = unnormalized object
#normalized = normalized object

testNormalizationGroups = function(unnormalized, normalized){
  coefV1 = sd(unnormalized$cellValsFrame$colSums)/mean(unnormalized$cellValsFrame$colSums)
  coefV2 = sd(normalized$colSums)/mean(normalized$colSums)
  return(1-(coefV2/coefV1))
}


