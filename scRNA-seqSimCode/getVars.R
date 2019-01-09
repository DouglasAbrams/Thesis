#calculates the variance in a set of clustering by calculating an anova table
#used in SNNcliq
#reduced_Counts = counts matrix post PCA
#clust = set of clusterings

getVariance = function(reduced_counts,clust){
  reduced_counts=cbind(reduced_counts,clust)
  
  anova1 = anova(lm(formula=reduced_counts[,1]~reduced_counts[,3],data=as.data.frame(reduced_counts)))
  anova2 = anova(lm(formula=reduced_counts[,2]~reduced_counts[,3],data=as.data.frame(reduced_counts)))
  
  p = anova1["reduced_counts[, 3]","Pr(>F)"] + anova2["reduced_counts[, 3]","Pr(>F)"]
  return(p)
}
