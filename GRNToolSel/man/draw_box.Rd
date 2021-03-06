% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_box.R
\name{draw_box}
\alias{draw_box}
\title{draw_box}
\usage{
draw_box(dat, name = "AUROC")
}
\arguments{
\item{dat}{The data.frame with the name of algorithm and the value of auc of aucpr.}

\item{name}{AUC or AUCPR}
}
\value{
The plot of AUC or AUCPR by different algorithm of bootstrap
}
\description{
Show the plot of AUC or AUCPR by different algorithms of bootstrap.
}
\examples{
set.seed(123)
A <- matrix(rnorm(200,0,1),20,10)
B <- matrix(rnorm(1000,0,1),100)
X <- matrix(0,100,20)
s <- c(1:10)
TF <- paste0('G',1:20)
gene <- paste0('G',1:100)
HGSmat <- matrix(0,100,20)
HGSmat[sample(2000,500)] <- 1
algorithm <- c('Pearson','ARACNE.a','PCIT','MRNETB','OMP')
conlist <- grn_main(A=A,B=B,X=X,seq = s,gene=gene,TF=TF,algorithm = algorithm,
                    HGSmat = HGSmat,file=FALSE,verbose=TRUE)
dat <- data.frame(Algorithm = algorithm,AUROC = conlist$AUROC)
p <- draw_box(dat,'AUROC')
print(p)
}
