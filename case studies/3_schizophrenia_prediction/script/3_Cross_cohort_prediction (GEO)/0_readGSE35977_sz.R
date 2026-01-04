rm(list=ls())

library(GEOquery)
library(affy)
library(hugene10stv1cdf)
library(hugene10stv1probe)
library(hugene10stprobeset.db)
library(hugene10sttranscriptcluster.db)


#unzip raw data
raw_data_directory = "../../raw_data/GEO/GSE35977_RAW/"
#cels <- list.files(raw_data_directory, pattern = "[gz]")
#head(cels)
#sapply(paste(raw_data_directory, cels, sep="/"), gunzip)


#read raw data
cels <- list.files(raw_data_directory, pattern = "CEL")
head(cels)
raw.data <- ReadAffy(verbose=TRUE, filenames=cels, celfile.path=raw_data_directory)


#calc expression level and normalization
expDat=exprs(rma(raw.data))
probes=rownames(expDat)
head(probes)

#check the anotation file name
raw.data
gene.symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
head(gene.symbols)
rownames(expDat) <- gene.symbols
expDat <- expDat[-which(is.na(gene.symbols)), ]
cels <- colnames(expDat)
head(cels)
cels <- sapply(cels, function (x) unlist(strsplit(x, "[.]"))[1])
head(cels)
colnames(expDat) <- cels




#get meta diagnosis data
gse <- getGEO("GSE35977", GSEMatrix=FALSE, destdir=raw_data_directory)
gsm <- GSMList(gse)
all(names(gsm)==colnames(expDat))
pheDat <- as.data.frame(names(gsm))
head(pheDat)


#add disease status
pheno <- unlist(lapply(gsm, function(x) x@header$characteristics_ch1[2]))
head(pheno)
pheno <- sapply(pheno, function(x) strsplit(x,split=" ")[[1]][3])
head(pheno)

pheno[pheno=="unaffected"] <- 0
pheno[pheno=="schizophrenia"] <- 1
pheno[pheno=="bipolar"] <- -1
pheno[pheno=="depression"] <- -2
pheno[pheno=="NA"] <- -100
head(pheno)
pheno <- as.numeric(pheno)
head(pheno)
pheDat <- cbind(pheDat,pheno)

#add age
ages <- sapply(gsm, function(x) x@header$characteristics_ch1[3])
head(ages)
ages <- sapply(ages, function(x) unlist(strsplit(x, " "))[2])
head(ages)
ages <- as.numeric(ages)
head(ages)
pheDat <- cbind(pheDat,ages)

#add gender
gender <- sapply(gsm, function(x) x@header$characteristics_ch1[5])
head(gender)
gender <- sapply(gender, function(x) unlist(strsplit(x, " "))[2])
head(gender)
gender[gender=='M'] <- 'Male'
gender[gender=='F'] <- 'Female'
head(gender)
pheDat <- cbind(pheDat,gender)

#add brain pH
ph <- sapply(gsm, function(x) x@header$characteristics_ch1[6])
head(ph)
ph <- sapply(ph, function(x) unlist(strsplit(x, " "))[2])
head(ph)
ph <- as.numeric(ph)
head(ph)
pheDat <- cbind(pheDat,ph)

#add post mortem interval
pmi <- sapply(gsm, function(x) x@header$characteristics_ch1[4])
head(pmi)
pmi <- sapply(pmi, function(x) unlist(strsplit(x, " "))[4])
head(ph)
pmi <- as.numeric(pmi)
head(ph)
pheDat <- cbind(pheDat,pmi)

#add rin
rin <- NA
pheDat <- cbind(pheDat,rin)
head(pheDat)

#Na Samples From X And Y
other <- as.character(pheDat[as.numeric(pheDat[,2])<0,1])
head(other)
expDat <- expDat[,-which(is.element(colnames(expDat), other))]
pheDat <- pheDat[-which(is.element(pheDat[,1], other)),]
dim(expDat)
dim(pheDat)




#save data
saveRDS(expDat, "../../preprocessed_data/GEO/expDat_sz_GSE35977.rds")
saveRDS(pheDat, "../../preprocessed_data/GEO/pheDat_sz_GSE35977.rds")
