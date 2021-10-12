Aclassify <- function(x_train,
                      x_test) {

    model <- keras::keras_model_sequential()
    model %>%
        keras::layer_dense(units = 256, activation = "relu", input_shape = ncol(x_train)) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_dense(units = 128, activation = "relu") %>%
        keras::layer_dense(units = 64, activation = "relu") %>%
        keras::layer_dense(units = 128, activation = "relu") %>%
        keras::layer_dense(units = 256, activation = "relu") %>%
        keras::layer_dense(units = ncol(x_train))
    #compile our model, using the mean squared error loss and the Adam optimizer for training
    model %>% keras::compile(
        loss = "mean_squared_error",
        optimizer = "adam",
        metrics = c('accuracy')
    )
    early_stopping <- callback_early_stopping(patience = 5)

    options(keras.view_metrics = FALSE)

    model %>% keras::fit(
        x = x_train,
        y = x_train,
        epochs = 100,
        batch_size = 32,
        validation_split = 0.2,
        callbacks = list(early_stopping)
    )
    #After training we can get the final loss for the test set by using the evaluate() fucntion.
    loss <- keras::evaluate(model, x = x_test, y = x_test)
    #Making predictions
    pred_test <- predict(model, x_test)
    mse_test <- apply((x_test - pred_test)^2, 1, sum)
    hc0 <- kmeans(mse_test, centers = 2)
    rec_err <- abs(x_test - pred_test)

    thiscluster <- hc0$cluster

    stopsign = FALSE
    failsign = FALSE
    niter = 1
    while(stopsign == FALSE & niter <= 10) {
        message(paste0("Feature selection: iteration", niter, "\n"))
        if(mean(mse_test[thiscluster == 1]) < mean(mse_test[thiscluster == 2])) {
            normalidx = 1
            tumoridx = 2
        } else {
            normalidx = 2
            tumoridx = 1
        }
        tumor_recErr <- rec_err[thiscluster == tumoridx, ]
        normal_recErr <- rec_err[thiscluster == normalidx, ]
        p=dim(tumor_recErr)[2]
        if(is.null(p)|length(p)==0){
            stopsign = TRUE
            failsign = TRUE
        }else if(is.na(p)){
            stopsign = TRUE
            failsign = TRUE
        }
        fullError <- rbind(tumor_recErr, normal_recErr)
        if(is.null(nrow(tumor_recErr)) | is.null(nrow(normal_recErr))) {
            message("Cannot detect any novel cells! Returning potential cluster based on the last iteration.")
            break
            stopsign = TRUE
            failsign = TRUE
            break
        } else {
            cttres <- genefilter::colttests(as.matrix(fullError), fac=factor(c(rep(1, nrow(tumor_recErr)), rep(0, nrow(normal_recErr)))))
            goodidx <- order(cttres$p.value)[1:500]

            rec_err2 <- rec_err[, goodidx]
            hc <- hclust(dist(rec_err2))
            memb <- cutree(hc, k = 2)

            message(paste0("#Cluster member change = ", sum(abs(thiscluster - memb)), "\n"))
            if(sum(abs(thiscluster - memb))<=2) {
                thiscluster <- memb
                stopsign = TRUE
            } else {
                thiscluster <- memb
            }
            niter = niter + 1
        }
    }

    if(!failsign & sum(thiscluster == 1)>1 & sum(thiscluster == 2)>1) {
        prof1 <- apply(x_train, 2, median)
        prof2_class1 <- apply(x_test[thiscluster == 1,], 2, median)
        prof2_class2 <- apply(x_test[thiscluster == 2,], 2, median)

        if(cor(prof1, prof2_class1) < cor(prof1, prof2_class2) | sum(colMeans(abs(rec_err[thiscluster ==1, ]))) > sum(colMeans(abs(rec_err[thiscluster ==2, ])))) {
            memb=thiscluster
            thiscluster[memb==1]=2
            thiscluster[memb==2]=1
        }

        finalcluster <- thiscluster - 1

        return(list(finalcluster = finalcluster,
                    rec_err = rec_err,
                    selectedGene = colnames(rec_err[, goodidx]),
                    failsign = failsign))
    } else {
        finalcluster <- thiscluster - 1

        return(list(finalcluster = finalcluster,
                    rec_err = rec_err,
                    failsign = failsign))
    }
}
