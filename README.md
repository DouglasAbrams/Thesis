# Thesis
Resolving UMI sequencing errors through machine learning and network-based approaches to improve the accuracy of rare cell type identification in single cell RNA-sequencing data.

### Abstract

  This study aims to develop new machine learning and network-based methods to  resolve UMI sequencing errors in single cell RNA-seq experiments. Until recently, UMI sequencing errors were poorly defined and addressed in ad hoc methods. Formalized methods to resolve them, introduced in Smith et. al. (2017), were successful, but imperfect. In helping to advance methods resolving UMI sequencing errors, we hope to improve the accuracy and effectiveness of single cell RNA-sequencing, which has and will continue to broaden our understanding of cellular taxonomy and heterogeneity.

### Background and Significance

  Unique Molecular Identifiers (UMIs), short oligonucleotides used to distinguish all of the cDNA molecules in a sequencing library, have recently become widely used in single cell RNA-sequencing (scRNA-seq). In such experiments, UMIs are used to reduce mRNA quantification errors due to PCR amplification biases by distinguishing correct molecular duplicates from erroneous ones. The use of UMIs in scRNA-seq is particularly crucial because of the high amount of amplification needed to make reads detectable during sequencing. Moreover, one of the primary goals of scRNA-seq is to identify rare and low-abundance cell types present in the sequencing library. This task is often made easier by UMIs as they can detect and reduce the numbers of abundant, over-amplified molecules, allowing for better detection of lower-abundance molecules, such as those from rare cell groups.


  Until recently, UMI sequencing errors, which have been shown to reduce the accuracy of downstream analyses, have been addressed by a variety of ad hoc methods, ranging from the removal of all molecules with extremely low-andundant UMIs to network based methods which seek to collapse networks of molecules containing nearly identical UMIs. In the first paper published directly addressing the question of resolving UMI sequencing errors, the authors proposed three different methods. The best method, which created and trimmed directional networks of reads with similar UMIs based on edit distance, improved the pearson correlation between ERCC spike-ins and read counts from 0.86 to 0.91 and reduced the number of miss-identified cells (identification was done through clustering) from 7 to 5. Though extremely successful, this method needs to be improved upon. While the authors were unable to distinguish the misidentified cells by type of error (read sequencing errors, clustering errors, limitations in resolving UMI sequencing errors), many scRNA-seq datasets have rare cell populations as low as 5 cells, and so it is important to further reduce rates of cell type misidentification.   
  
  
In this thesis, I plan to develop new machine learning and network-based methods to resolve UMI sequencing errors and to assess the performance of those that have already been proposed on both simulated scRNA-seq datasets and on real, publically available scRNA-seq datasets.
