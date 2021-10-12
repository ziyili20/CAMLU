fullAnno <- function(labels,
                     x_train,
                     y_train,
                     x_test,
                     method,
                     nfeature_svm = 1000) {
    if (length(dim(x_test[, labels == 0])) == 0) {
        newtest <- matrix(x_test[, labels == 0], nrow = length(x_test[, labels == 0]), ncol = 1)
    } else {
        newtest <- x_test[, labels == 0]
    }
    pred.out <- GetLabels_SVM(train_count = x_train,
                          train_label = y_train,
                          test_count = newtest,
                          nfeature_svm = nfeature_svm)
    fulllabel <- rep(NA, length(labels))
    fulllabel[labels == 0] <- pred.out
    fulllabel[labels == 1] <- "Unknown"

    return(fulllabel)
}
