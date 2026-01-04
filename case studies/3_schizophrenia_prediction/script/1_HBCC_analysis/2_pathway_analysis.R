
rm(list=ls())
gc()



fit = readRDS(file="fit_C=1_C2=1.rds")
backgroundGenes=rownames(fit$ws[[1]])
selectedGenes=backgroundGenes[rowMeans(fit$ws[[1]])!=0]

library(clusterProfiler)
library(org.Hs.eg.db)
pw <- enrichGO(gene=selectedGenes, universe=backgroundGenes,keyType="SYMBOL", OrgDb=org.Hs.eg.db, ont="ALL", pAdjustMethod = "fdr", 
               minGSSize=10, maxGSSize=200, pvalueCutoff=0.05, qvalueCutoff=0.05)
pw@result[,c("Description", "p.adjust")]