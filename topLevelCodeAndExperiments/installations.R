#installations

  
install.packages("devtools")
install.packages("bioconductor")
#import code from Github currently saved in a package
#devtools::install_github("DouglasAbrams/polyester")
#Biocundoctor
source("https://bioconductor.org/biocLite.R")
biocLite("splatter") #for scRNA-seq DE de-novo generation
biocLite("Biostrings") #non-installed dependency of polyester
install.packages("R.utils")?imp
#R.utils::sourceDirectory("/Volumes/Personal/dabrams/thesis/polyester-master/R") #path for colby macs
R.utils::sourceDirectory("/dabrams/thesis/polyester-master/R") #path for personal computer
r