rm(list=ls())
gc()

library(sva)
library(limma)
library(lumi)
library(affy)
library(annotate)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)

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

plotQC <- function(data, tag, path) {
  
  summary(data, 'QC')
  
  pdf(sprintf("%s/%s_norm_density.pdf",path,tag))
  plot(data, what='density')
  dev.off()
  
  pdf(sprintf("%s/%s_norm_boxplot.pdf",path,tag))
  plot(data, what='boxplot')
  dev.off()
  
  pdf(sprintf("%s/%s_norm_CDF.pdf",path,tag))
  plotCDF(data, reverse=TRUE)
  dev.off()
  
  pdf(sprintf("%s/%s_norm_pair.pdf",path,tag))
  pairs(data)
  dev.off()
  
  pdf(sprintf("%s/%s_norm_MA.pdf",path,tag))
  plot(data, what='MAplot')
  dev.off()
  
  pdf(sprintf("%s/%s_norm_cv.pdf",path,tag))
  plot(data, what='cv')
  dev.off()
  
  pdf(sprintf("%s/%s_norm_sample_relation.pdf",path,tag))
  plotSampleRelation(data, method='mds')
  dev.off()
}

detect_direction <- function(detection,expression){
  tmp <- expression[,1]
  # assumption is : a high detection value should have a low expression (intensity) value. 
  low <- tmp[!is.na(tmp)][which.max(detection[!is.na(tmp),1])]  
  # assumption is : a low detection value should have a high expression (intensity) value.
  high <- tmp[!is.na(tmp)][which.min(detection[!is.na(tmp),1])] 
  if (low > high) { # if assumption is not correct
    # detection <- 1 - detection # correct the detection data.table
    return(FALSE) # an expressed probe has a big detection value (close to 1). good quality expression : detectionTh > 0.9# 
  } 
  return(TRUE) # an expressed probe has a small detection value (close to 0). good quality expression : detectionTh < 0.0# 
  # see for more details : detectCall method of lumi : https://rdrr.io/bioc/lumi/src/R/detectionCall.R
}

#1, dlpfc isolation
#2, background correction
#3, quantile normalization
#4, log2 transformation
#5, remove low-expressed probes
#6, remove outliers
#7, sva


output_path <- "./"
qc_path <- "./QC_path"

raw <- readRDS("../../raw_data/HBCC/HumanHT-12_v4_lumi.rds")
meta <- read.table('../../raw_data/HBCC/metaData')

#filtering
metaExp <- meta[match(raw$sampleID, meta$SAMPLE_ID),]
metaExp <- metaExp[metaExp$HISTOLOGICAL_TYPE == "DLPFC",]
metaExp <- metaExp[metaExp$Diagnosis %in% c("Control"),] 	#extract control subjects
metaExp <- metaExp[metaExp$Race %in% c("African American","Caucasian"),]
metaExp$Diagnosis <- as.character(metaExp$Diagnosis)
metaExp=metaExp[!is.na(metaExp$AgeDeath), ]# remove 33
metaExp[metaExp$Smoker %in% "not filled in","Smoker"]="Unknown"

#missingness imputation
metaExp$pH[is.na(metaExp$pH)] <- mean(metaExp$pH[!is.na(metaExp$pH)])
metaExp$PMI[is.na(metaExp$PMI)] <- mean(metaExp$PMI[!is.na(metaExp$PMI)])
rownames(metaExp) <- metaExp$SAMPLE_ID

dim(raw)
dim(metaExp)
rawExp <- raw[, which(is.element(sampleNames(raw), metaExp$SAMPLE_ID))]
dim(rawExp)
dim(metaExp)
sum(metaExp$SAMPLE_ID!=sampleNames(rawExp))

# QC before Norm 
randomSamplesBeforeNorm <- rawExp[,sample(ncol(rawExp), 10)]
plotQC(randomSamplesBeforeNorm,"before",qc_path)

# normalization
lumi.N.Q <- lumiExpresso(rawExp, normalize.param = list(method='quantile'), 
                         varianceStabilize.param = list(method="log2", simpleOutput = FALSE), QC.evaluation=TRUE)

# filtering
expressionData=exprs(lumi.N.Q)
pcountA<-detectionCall(lumi.N.Q, Th=0.01)
hist(pcountA)
dim(expressionData)
expressionData<-expressionData[pcountA>ncol(lumi.N.Q)/2,]

lumi.N.Q <- lumi.N.Q[rownames(expressionData),]
dim(lumi.N.Q)

randomSamplesAfterNorm <- lumi.N.Q[,intersect(colnames(randomSamplesBeforeNorm),colnames(lumi.N.Q))]

#QC after Norm
plotQC(randomSamplesAfterNorm,"after",qc_path)

#annotation
lumi.N.Q <- addNuID2lumi(lumi.N.Q, annotationFile="../../raw_data/HBCC/HumanHT-12_V4_0_R2_15002873_B.txt")
symbol <- getSYMBOL(featureNames(lumi.N.Q), 'lumiHumanAll.db')
head(symbol)

expressionData <- exprs(lumi.N.Q)
gene <- as.vector(symbol[rownames(expressionData)])
expressionData <- cbind(gene,as.data.frame(expressionData))
expressionData <- expressionData[!is.na(expressionData[,1]),]
expressionData <- expressionData[order(apply(expressionData[,-1], 1, median),decreasing=T),]
expressionData <- expressionData[!duplicated(expressionData[1]), ]
rownames(expressionData) <- expressionData[,1]
expressionData$gene <- NULL

expressionData <- t(expressionData)  # rows: samples, columns: features

dim(expressionData)
outliers_removed <- remove_outliers(expressionData, metaExp)
metaExp <- outliers_removed$data_cov
data_1 <-outliers_removed$data
dim(data_1)

#sva analysis
mod = model.matrix(~AgeDeath+as.numeric(pH)+as.numeric(RIN)+as.numeric(PMI)+as.factor(Race) + as.factor(Sex), data=metaExp)
mod0 = model.matrix(~as.numeric(pH)+as.numeric(RIN)+as.numeric(PMI)+as.factor(Race)+as.factor(Sex)+1, data=metaExp)
ns.sv = 5
SVAsol = sva(dat=t(data_1), mod=mod, mod0=mod0, n.sv = ns.sv)
svs = SVAsol$sv

colnames(svs) <- paste0("sv", 1:5)
metaExp <- cbind(metaExp, svs)

# calculate PCA
ns.pca = 5
pctest = prcomp(data_1, scale=TRUE)
pcas <- pctest$x[,1:ns.pca]
metaExp <- cbind(metaExp, pcas)

# Look at PCA plot
pdf(sprintf("%s/before_res_PC.pdf",qc_path))
plot(pctest$x[,1], pctest$x[,2])
dev.off()

#residualize
covariates= metaExp[, c("pH", "RIN", "PMI", "Race", "Sex", paste0("sv", 1:5), paste0("PC", 1:5))]
all_res=apply(data_1, 2, function(x)lm(x~ ., data = covariates)$res)

colnames(metaExp)[which(names(metaExp) == "Diagnosis")] <- "diagnosis"
colnames(metaExp)[which(names(metaExp) == "Sex")] <- "sex"
colnames(metaExp)[which(names(metaExp) == "AgeDeath")] <- "age"
colnames(metaExp)[which(names(metaExp) == "Race")] <- "race"

# calculate PCA
ns.pca = 10
pctest = prcomp(all_res, scale=TRUE)
pcas <- pctest$x[,1:ns.pca]

# Look at PCA plot
pdf(sprintf("%s/after_res_PC.pdf",qc_path))
plot(pctest$x[,1], pctest$x[,2])
dev.off()


saveRDS(all_res, file="../../preprocessed_data/HBCC/gene_expression_predict_age.rds")
saveRDS(metaExp[,c("SAMPLE_ID", "sex", "age")], "../../preprocessed_data/HBCC/y_matrix_predict_age.rds" )
