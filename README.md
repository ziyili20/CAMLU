# Identify novel cells using CAMLU

`CAMLU` is an R package that provides an autoencoder based method for annotating cell types from scRNA-seq data. The function can identify unknown cells with the input training data. It also can annotate the full lists of cell types with consideration of unknown cell types. This vignette introduces the CAMLU function and the things it can do for you. CAMLU was created by Ziyi Li, Yizhuo Wang, Irene Ganan-Gomez, Simona Colla and Kim-Anh Do, and is now maintained by Ziyi Li (zli16@mdanderson.org).

<img width="1070" alt="Sche" src="https://user-images.githubusercontent.com/24440640/161843184-cb90a9a9-1264-408d-aa98-80c849af0082.png">

## Installation


To install this package, start R (version "3.6" or higher) and enter:

```{r install, message=FALSE, warning=FALSE}
install.packages("devtools")
library(devtools)

install_github("ziyili20/CAMLU", build_vignettes = TRUE)
library(CAMLU)
```

Note that CAMLU relies on the R keras package (https://cran.r-project.org/web/packages/keras/index.html). You can install it through:
```
install.packages('keras')
```

### How to get help for CAMLU

Any questions should be posted
to the GitHub Issue section of CAMLU
homepage at https://github.com/ziyili20/CAMLU/issues.

### CAMLU function

The main function has 6 inputs:

-   x_train: a scRNA-seq dataset as the training data.
-   x_test: another scRNA-seq dataset as the testing data, note that the function identifies the cell type that is not present in the training data, only exists in the testing data.
-   full_annotation: an optional step to annotate the full range of cell types; the default is FALSE.
-   y\_train: full labels of the training data; should be provided if full_annotation is set up to TRUE.
-   ngene: the number of top variable genes in the feature selection step for identifying unknown cells; if missing, it will be defaulted to 5000.
-   lognormalize: the log10 normalization step applying to training and testing datasets; the default is TRUE unless specify.


```{r run1, eval = TRUE, message = FALSE}
## Load data
library(CAMLU)
data(PBMC_tumor_simulation_data)
## Distinguish novel cells from known cell types
label_01 <- CAMLU(x_train = simdata$fdata_train,
                  x_test = simdata$fdata_test,
                  ngene=5000,lognormalize=TRUE)
## Display the classification table
truelabel <- simdata$labels_2
print(table(label_01, truelabel))
```

In the results, label\_01 = 0 represent the labels for the normal cells and the label\_01 = 1 are the tumor cells. Let's visualize the cells in the TSNE plot.

```{r start, message = TRUE}
library(SingleCellExperiment)
library(scater)
sce <- SingleCellExperiment(list(counts=simdata$fdata_test),
                            colData=DataFrame(celltype=as.factor(truelabel),
predicted = as.factor(label_01))
                            )
sce <- logNormCounts(sce) ## add log normalized counts
sce <- runTSNE(sce, perplexity=10)
par(mfrow = c(1,2))
plotTSNE(sce, colour_by = "celltype")
plotTSNE(sce, colour_by = "predicted")
```

In addition to annotating the normal/tumor cells, we can also annotate the full lists of cell types by specifying `full_annotation = TRUE`. 

```{r full, message=TRUE, fig.height=5, fig.width=7}
## Enable full annotation function: assign the full spectrum of labels
label_all <- CAMLU(x_train = simdata$fdata_train,
                  x_test = simdata$fdata_test,
                  y_train = simdata$labels_train,
                  full_annotation = TRUE,
                  ngene=3000,
                  lognormalize=TRUE)
## Display the classification table
truelabel_all <- simdata$labels_test
truelabel_all[truelabel_all == "HNCC"] <- "Unknown"
table(label_all$label_full, truelabel_all)
## Calculate the accuracy
sum(label_all$label_full == truelabel_all)/ length(truelabel_all)
```
Let's visualize the annotation results again.

```{r vis2, message = TRUE}
sce <- SingleCellExperiment(list(counts=simdata$fdata_test),
                            colData=DataFrame(celltype=as.factor(simdata$labels_test),
                                              predicted = as.factor(label_all$label_full))
)
sce <- logNormCounts(sce) ## add log normalized counts
sce <- runTSNE(sce, perplexity=10)
par(mfrow = c(1,2))
plotTSNE(sce, colour_by = "celltype")
plotTSNE(sce, colour_by = "predicted")

```
## Reference

Zheng, G. X., Terry, J. M., Belgrader, P., Ryvkin, P., Bent, Z.W.,Wilson, R., Ziraldo, S. B., Wheeler, T. D., McDermott, G. P., Zhu, J., *et al.* (2017). Massively parallel digital transcriptional profiling of single cells. *Nature communications*, 8(1), 1--12.

