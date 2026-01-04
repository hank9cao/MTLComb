rm(list=ls())
gc()

library(Biobase)
library(GenomicRanges)
library(DESeq2)
library(rafalib)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sva)
library(dplyr)
library(edgeR)
library(EnsDb.Hsapiens.v86)
library(scales)
library(affy)
library(remotes)
library(RNAseqQC)
library(ensembldb)
library(ggplot2)
library(purrr)
library(tidyr)
library(tibble)
library(magrittr)

#remove outliers
remove_outliers <- function(data_pca, data_pca_cov) {   # data_pca rows : samples, columns: features
  curPCA = prcomp(data_pca, scale=TRUE)$x[,1:2]
  out1 = which( abs((curPCA[,1] - mean(curPCA[,1]))/sd(curPCA[,1])) > 4 )
  out2 = which( abs((curPCA[,2] - mean(curPCA[,2]))/sd(curPCA[,2])) > 4 )
  outlier = unique(sort(c(out1, out2)))
  if(length(outlier)>0){
    data_new  <- data_pca[-outlier, ]
    data_cov  <- data_pca_cov[-outlier, ]
  } else {
    data_new <- data_pca
    data_cov <- data_pca_cov
  }
  result = list()
  result$data = data_new
  result$data_cov = data_cov
  return(result) 	# output rows: samples, columns: features
}

plotQC <- function(data,path,suffix,threshold,samples) {
  
  myColors <- hue_pal()(4)
  pdf(sprintf("%s/boxplot_%s.pdf",path,suffix))
  boxplot(log2(data[,samples]+0.1), xlab="", ylab="Log2 RPKM")
  dev.off()
  
  
  pdf(sprintf("%s/density_%s.pdf",path,suffix))
  plotDensity(log2(data+0.1), col=rep(myColors, each=3),
              lty=c(1:ncol(data)), xlab='Log2(rpkm)',
              main='Expression Distribution')
  dev.off()
  
  pdf(sprintf("%s/hist_%s.pdf",path,suffix))
  hist(rowSums(data > threshold), breaks = 50)
  dev.off()
}

qc_path <- "./QC_path/"
load("../../raw_data/Brainseq/rse_gene_unfiltered.Rdata")

# filtering
rse_gene <- rse_gene[, rse_gene$Age >= 18]
rse_gene <- rse_gene[, rse_gene$Age <= 65]
rse_gene <- rse_gene[, rse_gene$Region == "DLPFC"]
rse_gene <- rse_gene[, is.element(rse_gene$Race, c("AA", "CAUC"))]
dim(rse_gene) #320

# convert gene symbols to gene names
ensemblIds <- unlist(lapply(stringr::str_split(rownames(rse_gene), "[.]"), "[[", 1))
geneSymbols <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensemblIds, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(rse_gene) <- geneSymbols$SYMBOL[match(ensemblIds,geneSymbols$GENEID)]


#Aggregate genes
rpkm<-assays(rse_gene)$rpkm
rpkm_aggregated <- rpkm[rowSums(rpkm) > 0,]
rpkm_aggregated <- rpkm_aggregated[!duplicated(rownames(rpkm_aggregated)), ] #<- difference? Nor keep one with maxed median

#replace label 
meta <- colData(rse_gene)
label = replace(meta$Dx, meta$Dx == 'Control', -1)
label = replace(label, label == 'Schizo', 1)
meta$Dx = as.integer(label)



stopifnot(identical(colnames(rpkm_aggregated),rownames(meta)))

######### only demographic visualisation ############ 
case_group <- meta[meta$Dx == 1,]
control_group <- meta[meta$Dx == -1,]

pdf(sprintf("%s/case_control.pdf",qc_path))
slices <- c(dim(case_group)[1], dim(control_group)[1])
lbls <- c(paste0("case ",dim(case_group)[1]), paste0("control ",dim(control_group)[1]))
pie(slices, labels = lbls, main="Case/Control")
dev.off()

pdf(sprintf("%s/male_female_case.pdf",qc_path))
slices <- c(sum(case_group$Sex == "F"), sum(case_group$Sex == "M"))
lbls <- c(paste0('Female ',sum(case_group$Sex == "F")), paste0('Male ',sum(case_group$Sex == "M")))
pie(slices, labels = lbls, main="Case")
dev.off()

pdf(sprintf("%s/male_female_control.pdf",qc_path))
slices <- c(sum(control_group$Sex == "F"), sum(control_group$Sex == "M"))
lbls <- c(paste0('Female ',sum(control_group$Sex == "F")), paste0('Male ',sum(control_group$Sex == "M")))
pie(slices, labels = lbls, main="Control")
dev.off()

pdf(sprintf("%s/age_dist_case.pdf",qc_path))
hist(case_group$Age)
dev.off()
pdf(sprintf("%s/age_dist_control.pdf",qc_path))
hist(control_group$Age)
dev.off()
#####################################################

threshold = 1

randomSamples <-sample(1:ncol(rpkm_aggregated),10)

plotQC(rpkm_aggregated, qc_path, "before", threshold, randomSamples)

normalized_counts <- rpkm_aggregated[rowSums(rpkm_aggregated > threshold)> 75,] #? <- how the 75 come out

plotQC(normalized_counts, qc_path, "after", threshold, randomSamples)

bs_pca = remove_outliers(t(normalized_counts), meta)

data = bs_pca$data		# rows: samples, columns: features
#317 18815
meta = bs_pca$data_cov

age2 <- meta$Age^2  
meta <- cbind(meta, age2)

stopifnot(identical(rownames(meta), rownames(data)))

#cell proportion estimation
load("../../raw_data/Brainseq/methprop_pd.Rdata")
cells = c("Astrocytes","Oligodendrocytes","Microglia","Endothelial","Neurons","OPC")
meta = cbind(meta,pd[rownames(meta),cells])

#run sva
covariates= c("Age", "Race", "RIN", "age2", "Sex", cells)
str(meta[,c("Dx", covariates)])
mod = model.matrix(~as.factor(Dx) + Age + Race + as.numeric(RIN) + age2 + Sex + Astrocytes + Oligodendrocytes + Microglia + Endothelial + Neurons + OPC, 
                   data = meta)
mod0 = model.matrix(~1 + Age + Race + as.numeric(RIN) + age2 + Sex + Astrocytes + Oligodendrocytes + Microglia + Endothelial + Neurons + OPC, 
                    data = meta) #? why not cells
n.sv = 5
SVAsol = sva(dat = t(data), mod = mod, mod0 = mod0, n.sv = n.sv)
svs = SVAsol$sv
colnames(svs)=paste0("sv", 1:5)
meta = cbind(meta,svs)

#calculate PCAs
n.pca = 5
pctest = prcomp(data, scale=TRUE)
pcas <- pctest$x[,1:n.pca]
meta <- cbind(meta, pcas)

#Look at PCA before Residualization
pdf(sprintf("%s/before_res_PC.pdf",qc_path))
plot(pctest$x[,1], pctest$x[,2])
dev.off()

#residualize data
covariates <- meta[, c("Age","Race", "RIN","Sex", "age2", paste0("sv", 1:5), paste0("PC", 1:5),cells)] #<- no pc?
covariates$RIN=as.numeric(covariates$RIN)
corr_bs <- apply(data, 2, function(x){
  lm.obj <- lm(x~ ., data = covariates)$res
  #pmax(mean(x) + x - predict(lm.obj), 0) 
  #? <- This only preserv the positive residual, why?
  return(lm.obj)
})

#Look at PCA before Residualization
pctest = prcomp(corr_bs, scale=TRUE)
pcas <- pctest$x[,1:n.pca]

pdf(sprintf("%s/after_res_PC.pdf",qc_path))
plot(pctest$x[,1], pctest$x[,2])
dev.off()

#save
colnames(meta)[which(names(meta) == "Dx")] <- "diagnosis"
colnames(meta)[which(names(meta) == "Sex")] <- "sex"
colnames(meta)[which(names(meta) == "Age")] <- "age"
colnames(meta)[which(names(meta) == "Race")] <- "race"

saveRDS(corr_bs, file="../../preprocessed_data/Brainseq/gene_expression_predict_Dx.rds")
saveRDS(as.data.frame(meta[,c("sex", "age", "diagnosis")]), "../../preprocessed_data/Brainseq/y_matrix_predict_Dx.rds" )

