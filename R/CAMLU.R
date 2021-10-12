CAMLU <- function(x_train,
                  x_test,
                  full_annotation = FALSE,
                  y_train = NULL,
                  ngene = 5000,
                  lognormalize = TRUE) {

    pout <- PrepFunc(x_train = x_train,
                     x_test = x_test,
                     ngene = ngene,
                     lognormalize = lognormalize)
    aout <- Aclassify(x_train = pout$x_train,
                      x_test = pout$x_test)
    ## this is the label with known/unknown information
    label_known_unknown <- aout$finalcluster

    if (full_annotation) {
        if (is.null(y_train)) {
            stop("For full annotation function, please provide full labels for training data through y_train.")
        } else {
            fulllabel <- fullAnno(aout$finalcluster,
                                  x_train = x_train,
                                  y_train = y_train,
                                  x_test = x_test,
                                  method = "SVM")
            return(list(label_known_unknown = label_known_unknown,
                        label_full = fulllabel))
        }
    } else {
        return(label_known_unknown)
    }
}
