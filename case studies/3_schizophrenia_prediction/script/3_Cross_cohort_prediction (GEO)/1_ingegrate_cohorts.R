rm(list=ls())

sets <- c("GSE53987","GSE21138", "GSE35977")

X <- list()
Y <- list()

for (i in sets){
    print(i)
    #read x and y from each dataset
    tempX <- readRDS(paste("../../preprocessed_data/GEO/expDat_sz_", i, ".rds",sep=""))
    tempY <- readRDS(paste("../../preprocessed_data/GEO/pheDat_sz_", i, ".rds",sep=""))

    genes <- rownames(tempX)
    genes.uni <- unique(genes)

    #for x: average expression value for the gene
    X[[i]] <- sapply(genes.uni, function (x) {colMeans(tempX[is.element(genes,x), ,drop=FALSE])})

    #for x: show info
    print(length(genes))
    print(length(genes.uni))
    print(dim(X[[i]]))
    print(head(colnames(X[[i]])))
    print(head(rownames(X[[i]])))

    #for y: remove first column and attach rownames
    Y[[i]] <- tempY

    #double check the alignment of subject
    n1 <- rownames(X[[i]])
    n2 <- rownames(Y[[i]])
    print(all(n1==n2))

    print(dim(tempY))
    print(head(rownames(tempY)))
}

T=length(X)

#calc common genes
genes.uni <- colnames(X[[1]])
for (i in (2:T)){
  genes.uni <- intersect(genes.uni, colnames(X[[i]]))
}

for (i in (1:T)){
  X[[i]] <- X[[i]][,match(genes.uni, colnames(X[[i]]))]
  print(head(colnames(X[[i]])))
}


save(file='../../preprocessed_data/GEO/unnormalized_X_Y.rda', X, Y)



