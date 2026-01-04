rm(list=ls())
gc()



fits = readRDS(file="fit1_fit2_test.rds")
plot(rowMeans(fits[[1]]$ws[[1]]), rowMeans(fits[[2]]$ws[[1]]))
fit=fits[[1]]
backgroundGenes=rownames(fit$ws[[1]])
selectedGenes=names(which(rowMeans(fit$ws[[1]])!=0))
selectedGenes=names(sort(abs(rowMeans(fit$ws[[1]])), decreasing = T)[1:100])

library(clusterProfiler)
library(org.Hs.eg.db)
pw <- enrichGO(gene=selectedGenes, universe=backgroundGenes,keyType="SYMBOL", OrgDb=org.Hs.eg.db, ont="ALL", pAdjustMethod = "fdr", 
               minGSSize=10, maxGSSize=200, pvalueCutoff=0.05, qvalueCutoff=0.05)
pw@result[,c("Description", "p.adjust")]


#GO:0005903                            brush border 0.005935157
#GO:0098862 cluster of actin-based cell projections 0.007792799
