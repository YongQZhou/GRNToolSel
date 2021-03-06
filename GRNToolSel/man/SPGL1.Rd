% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPGL1.R
\name{SPGL1}
\alias{SPGL1}
\title{SPGL1 - spectral projected gradient method to solve L1 regularization problem}
\usage{
SPGL1(A, b, x0, sparsity, MaxIter = 200, epsilon = 1e-06)
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
This function is the process of SPGL1 algorithm aims to solve \eqn{l_1} regularization optimization model.
}
\details{
Spectral projected gradient method (SGPL1) was proposed by van den Berg and Friedlander (2008) to solve the \eqn{l_1} optimization problem:
\deqn{{\min}_{x \in R^n} ||Ax - b|{|_2^2} + \lambda ||x|{|_1}}
SPGL1 is indeed the projected gradient method with Barzilai-Borwein stepsize rule Barzilai and Borwein (1988), that is, first takes the gradient descent step with Barzilai-Borwein line-search scheme and then proceeds the projection onto the \eqn{l_1} norm ball at each iteration.
}
\references{
van den Berg, E. and Friedlander, M. P. (2008). "Probing the pareto frontier for basis pursuit solutions", SIAM Journal on Scientific Computing, 31(2), 890-912.
}
\author{
Yaohua Hu <mayhhu@szu.edu.cn> 

Xinlin Hu <ttxinlinhu@qq.com>
}
