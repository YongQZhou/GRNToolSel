% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ITAHard.R
\name{ITAHard}
\alias{ITAHard}
\title{ITAHard - iterative hard thresholding algorithm}
\usage{
ITAHard(A, b, x0, sparsity, MaxIter = 200, epsilon = 1e-06)
}
\arguments{
\item{A}{Gene expression data of transcriptome factors (i.e. basis function in machine learning).}

\item{b}{Gene expression data of target genes (i.e. observation in machine learning).}

\item{x0}{Gene expression data of Chromatin immunoprecipitation or zero vector (i.e. initial iterative point in machine learning).}

\item{sparsity}{Sparsity level of solution.}

\item{MaxIter}{Maximum iteration used in calculation; we set default as 200.}

\item{epsilon}{Stopping rule in algorithm where the rule is \eqn{||Ax-b|{|_2^2} < \epsilon}; we set default as 1e-6.}
}
\description{
This function is the process of ITA-Hard algorithm aims to solve \eqn{l_0} regularization optimization model, and it is an extension of ISTA.
}
\details{
Iterative hard thresholding algorithm (ITA-Hard) was introduced by Blumensath and Davies (2008), which is an extension of the ISTA for solving the \eqn{l_0} regularization problem:
\deqn{{\min}_{x \in R^n} ||Ax - b|{|_2^2} + \lambda ||x|{|_0}}
Similar to ISTA, the idea of ITA-Hard is to sequentially proceed a gradient descent step and a hard thresholding operator at each iteration.
}
\references{
Blumensath, T. and Davies, M. E. (2008). "Iterative thresholding for sparse approximations", Journal of Fourier Analysis and Applications, 14, 629-654.
}
\author{
Yaohua Hu <mayhhu@szu.edu.cn> 

Xinlin Hu <ttxinlinhu@qq.com>
}
