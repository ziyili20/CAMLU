PrepFunc <- function(x_train,
                     x_test,
                     ngene=5000,
                     lognormalize=TRUE) {

    ## filter genes a bit
    ## only keep genes that have row summation more than 10
    x_train<-x_train[rowSums(x_train>0)>=10,]
    x_test<-x_test[rowSums(x_test>0)>=10,]

    ## find overlapped genes between the training and testing datasets
    gene_name=intersect(rownames(x_train),rownames(x_test))
    x_test=x_test[gene_name,]
    x_train=x_train[gene_name,]

    message(paste0("Number of overlapped genes in training and testing dataset:", length(gene_name), "."))

    ## log-normalize data
    sce.data <- SingleCellExperiment::SingleCellExperiment(list(counts=x_train))
    sce.data <- scuttle::logNormCounts(sce.data)

    ## select highly variable genes
    dec.data <- scran::modelGeneVar(sce.data)
    data.var <- scran::getTopHVGs(dec.data, n = ngene)
    if(length(data.var)<ngene){
        message(paste0("Number of selected high variance gene is less than specified ngene:", length(data.var), "."))
    }

    x_train=x_train[data.var,]
    sce.data <- SingleCellExperiment::SingleCellExperiment(list(counts=x_train))
    sce.data <- scuttle::logNormCounts(sce.data)
    if(lognormalize==FALSE){
        x_train=t(sce.data@assays@data$counts)
    }else{
        x_train=t(sce.data@assays@data$logcounts)
    }

    x_test=x_test[data.var,]
    sce.data <- SingleCellExperiment::SingleCellExperiment(list(counts=x_test))
    sce.data <- scuttle::logNormCounts(sce.data)
    if(lognormalize==FALSE){
        x_test=t(sce.data@assays@data$counts)
    }else{
        x_test=t(sce.data@assays@data$logcounts)
    }

    return(list(x_train = x_train,
                x_test = x_test))
}
