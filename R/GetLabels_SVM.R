GetLabels_SVM <- function(train_count = NULL,
                          train_label = NULL,
                          test_count,
                          nfeature_svm = 1000){

    tsd <- apply(train_count, 1, sd)
    ooidx <- order(tsd, decreasing = TRUE)[1:nfeature_svm]

    train_count2 <- train_count[ooidx,]
    test_count2 <- test_count[ooidx,]

    uniqID <- unique(train_label)
    matlab <- match(train_label, uniqID)

    ## train
    svmout <- e1071::svm(t(as.matrix(train_count2)), as.factor(matlab), kernel = "linear")

    ## predict
    predlab <- predict(svmout, t(as.matrix(test_count2)))

    outlab <- uniqID[predlab]
    return(outlab)
}
