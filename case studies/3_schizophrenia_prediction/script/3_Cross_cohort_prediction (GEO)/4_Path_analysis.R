rm(list=ls())
gc()



fits = readRDS(file="fit1_fit2_test.rds")
plot(rowMeans(fits[[1]]$ws[[1]]), rowMeans(fits[[2]]$ws[[1]]))
fit=fits[[2]]
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





rm(list=ls())
gc()

fits = readRDS(file="fit_cvFit_all.rds")
fit=fits[[1]]
backgroundGenes=rownames(fit$ws[[1]])
selectedGenes=names(which(rowMeans(fit$ws[[1]])!=0))
selectedGenes=names(sort(abs(rowMeans(fit$ws[[1]])), decreasing = T)[1:500])

library(clusterProfiler)
library(org.Hs.eg.db)
pw <- enrichGO(gene=selectedGenes, universe=backgroundGenes,keyType="SYMBOL", OrgDb=org.Hs.eg.db, ont="ALL", pAdjustMethod = "fdr", 
               minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, qvalueCutoff=0.05)
pw@result[,c("Description", "p.adjust")]


# GO:0007268                 chemical synaptic transmission 0.025878397
# GO:0098916           anterograde trans-synaptic signaling 0.025878397
# GO:0099537                       trans-synaptic signaling 0.025878397
# GO:0009887                     animal organ morphogenesis 0.025878397
# GO:0030001                            metal ion transport 0.025878397
# GO:0006469 negative regulation of protein kinase activity 0.025878397
# GO:0007423                      sensory organ development 0.032789873
# GO:0099536                             synaptic signaling 0.033531818
# GO:0033673         negative regulation of kinase activity 0.046850028
# GO:0005244  voltage-gated monoatomic ion channel activity 0.003115691
# GO:0022832                 voltage-gated channel activity 0.003115691
# GO:0022836                         gated channel activity 0.007911141
# GO:0022839          monoatomic ion gated channel activity 0.007911141

