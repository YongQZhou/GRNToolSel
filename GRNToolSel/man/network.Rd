% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{network}
\alias{network}
\title{Visualization of network}
\usage{
network(grn, num = NULL, TF, gene, file = NULL)
}
\arguments{
\item{grn}{A matrix which is the weighted adjacency matrix of the inferred network by this algorithm or the result from the function of out_grn.}

\item{num}{The number of gene regulator imformation. Default: NULL.}

\item{TF}{The Transcription Factors (TFs) The default value NULL means that all the genes are used as candidate regulators. Default: NULL.}

\item{gene}{A vector of characters containing the target genes.}

\item{file}{A connection, or a character string naming the file to write to. The default not save. Default: NULL.}
}
\value{
A D3 JavaScript force directed network graph.
}
\description{
Show the gene regulator network which you provided. Visualization of gene regulatory network
}
\details{
A labelled undirected network is plotted.
}
\examples{
B <- matrix(rnorm(1000,0,1),100)
TF <- paste0('G',1:20)
gene <- paste0('G',1:100)
grn <- demo_other(B,'ARACNE.a',TF,gene)$grn
net <- network(grn,TF=TF,gene=gene)
}
