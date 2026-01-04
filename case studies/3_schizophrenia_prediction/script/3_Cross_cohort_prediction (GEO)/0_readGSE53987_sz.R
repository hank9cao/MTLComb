rm(list=ls())


library(GEOquery)
library(affy)
library(hgu133plus2.db)
library(annotate)


#unzip raw data
raw_data_directory = "../../raw_data/GEO/GSE53987_RAW/"
#cels <- list.files(raw_data_directory, pattern = "[gz]")
#head(cels)
#sapply(paste(raw_data_directory, cels, sep="/"), gunzip)


#read raw data
cels <- list.files(raw_data_directory, pattern = "CEL")
head(cels)
raw.data <- ReadAffy(verbose=TRUE, filenames=cels, celfile.path=raw_data_directory)


#calc expression level and normalization
expDat=exprs(rma(raw.data))
probes=row.names(expDat)
head(probes)

#check the anotation file name
gene.symbols <- getSYMBOL(probes,"hgu133plus2")
head(gene.symbols)
rownames(expDat) <- gene.symbols
expDat <- expDat[-which(is.na(gene.symbols)), ]
cels <- colnames(expDat)
head(cels)
cels <- sapply(cels, function (x) unlist(strsplit(x, "[.]"))[1])
head(cels)
cels <- sapply(cels, function (x) unlist(strsplit(x, "[_]"))[1])
head(cels)
colnames(expDat) <- cels


#get meta diagnosis data
gse <- getGEO("GSE53987", GSEMatrix=FALSE, destdir=raw_data_directory)
gsm <-GSMList(gse)
all(names(gsm)==colnames(expDat))
pheDat <- as.data.frame(names(gsm))
head(pheDat)

#add disease status
pheno <- unlist(lapply(gsm, function(x) x@header$title))
head(pheno)
pheno1 <- sapply(pheno, function(x) unlist(strsplit(x, "[_]"))[c(1,2)])
pheno <- rep(-99,ncol(pheno1))
head(pheno)
head(pheno1)
pheno[apply(pheno1,2,function(x){x[1]=="schiz" & x[2]=="pfc"})] <- 1
pheno[apply(pheno1,2,function(x){x[1]=="control" & x[2]=="pfc"})] <- 0
pheno <- as.numeric(pheno)
head(pheno)
pheDat <- cbind(pheDat,pheno)



#add age
ages <- sapply(gsm, function(x) x@header$characteristics_ch1[1])
head(ages)
ages <- sapply(ages, function(x) unlist(strsplit(x, " "))[2])
head(ages)
ages <- as.numeric(ages)
head(ages)
pheDat <- cbind(pheDat,ages)

#add gender
gender <- sapply(gsm, function(x) x@header$characteristics_ch1[2])
head(gender)
gender <- sapply(gender, function(x) unlist(strsplit(x, " "))[2])
head(gender)
gender[gender=='M'] <- 'Male'
gender[gender=='F'] <- 'Female'
pheDat <- cbind(pheDat,gender)

#add brain pH
ph <- sapply(gsm, function(x) x@header$characteristics_ch1[5])
head(ph)
ph <- sapply(ph, function(x) unlist(strsplit(x, " "))[2])
head(ph)
ph <- as.numeric(ph)
head(ph)
pheDat <- cbind(pheDat,ph)

#add post mortem interval
pmi <- sapply(gsm, function(x) x@header$characteristics_ch1[4])
head(pmi)
pmi <- sapply(pmi, function(x) unlist(strsplit(x, " "))[2])
head(pmi)
pheDat <- cbind(pheDat,pmi)

#add rin
rin <- sapply(gsm, function(x) x@header$characteristics_ch1[6])
head(rin)
rin <- sapply(rin, function(x) unlist(strsplit(x, " "))[2])
head(rin)
pheDat <- cbind(pheDat,rin)



#remove bipolar sample from X and Y
other <- as.character(pheDat[pheDat[,2]<0,1])
expDat <- expDat[,-which(is.element(colnames(expDat), other))]
pheDat <- pheDat[-which(is.element(pheDat[,1], other)),]
head(pheDat)


#save data
dim(expDat)
dim(pheDat)
saveRDS(expDat, "../../preprocessed_data/GEO/expDat_sz_GSE53987.rds")
saveRDS(pheDat, "../../preprocessed_data/GEO/pheDat_sz_GSE53987.rds")
