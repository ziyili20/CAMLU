\name{Cell annotation with the presence of novel cells}
\alias{CAMLU}
\title{CAMLU: Cell annotation with the presence of novel cells}
\usage{
CAMLU(x_train, x_test, full_annotation = FALSE, y_train = NULL, ngene = 5000, lognormalize = TRUE)
}
\description{
Annotate cells from single cell RNA-seq data with consideration of potential novel cells.
The function can either identify known/unknown cells, or annotate the full set of cell types.
}
\arguments{
\item{x_train}{A p1 by N1 matrix including scRNA-seq expression counts for training. Raw (unnormalized) counts should be provided. Rows are genes and columns are cells.}
\item{x_test}{A p2 by N2 matrix including scRNA-seq expression counts for annotation prediction. Raw (unnormalized) counts should be provided. Rows are genes and columns are cells.}
\item{full_annotation}{Whether to annotate the full lists of cell types or not. Default is FALSE, CAMLU only identify unknown cells.
If TRUE, the full vector of cell labels need to be provided through y_train. }
\item{y_train}{A N1-length vector including the cell labels of the N1 cells. If full_annotation is FALSE, y_train will not be used.}
\item{ngene}{Number of highly variable genes selected in the first analysis step. Default value is 5000.}
\item{lognormalize}{Whether the data should be log-normalized or not. Default is TRUE}
}
\examples{
### a toy example
x_train <- matrix(rnbinom(2000*1000, size = sample(1:50, 2000*1000, replace = TRUE), prob = 0.2), 2000, 1000)
x_test <- matrix(rnbinom(2000*1000, size = sample(1:50, 2000*1000, replace = TRUE), prob = 0.2), 2000, 1000)

rownames(x_train) <- paste0("gene", 1:2000)
rownames(x_test) <- paste0("gene", 1:2000)

label_known_unknown = CAMLU(x_train,
                            x_test,
                            full_annotation = FALSE,
                            y_train = NULL,
                            ngene = 5000,
                            lognormalize = TRUE)

y_train <- sample(paste0("Celltype", 1:8), ncol(x_train), replace = TRUE)
label_known_unknown = CAMLU(x_train,
                            x_test,
                            full_annotation = TRUE,
                            y_train = y_train,
                            ngene = 5000,
                            lognormalize = TRUE)
}
