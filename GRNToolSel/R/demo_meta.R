#' Demo meta - The demo of meta method
#' @description This is the main function to call meta algorithms for network inference
#'
#' @param algorithm The algorithm of inferring grn. We provide 25 algorithms in this function including
#' 'Pearson','Spearman','Kendall','WGCNA','BC3NET','C3NET','ARACNE.a','ARACNE.m','CLR','MRNET','MRNETB',
#' 'GeneNet','GENIE3.ET','GENIE3.RF','ENNET,'PCIT','ADMM', 'ADMMHalf',ADMMHard', 'ITAHalf', 'ITAHard',
#' 'ISTA', 'SPGL0', 'SPGL1', 'CoSaMP', 'OMP', 'FoBa' and 'LARS'. Defaule: c('Pearson','ARACNE.a','CLR')
#' @param TF The Transcription Factors (TFs) The default value NULL means that all the genes are used as candidate regulators. Default: NULL.
#' @param gene A vector of characters containing the target genes.
#' @param cl.cores The number of cores in computer that you want to utilize. Default: 2.
#' @param file A connection, or a character string naming the file to write to. The default not save the grn matrix. Default: NULL.
#' @param verbose If set to TRUE, a feedback on the progress of the progress of the calculations is given. Default: TRUE.
#' @param A Gene expression data of transcriptome factors (i.e. basis function in machine learning).
#' The class of A are required to be 'matrix' and the dimension of matrix A is m * n.
#' @param B Gene expression data of target genes (i.e. observation in machine learning).
#' The class of B are required to be 'matrix' and the dimension of matrix B is u * n.
#' @param X Gene expression data of Chromatin immunoprecipitation or zero matrix
#' (i.e. initial iterative point in machine learning).The class of X are required to be
#' 'matrix' and the dimension of matrix X is u * m.
#' @param s Sparsity level of solution. User can input a sequence of sparsity, i.e. 's <- c(1,2,3,4,5)'.
#' @param p If set to FALSE, this demo aims to use the algorithms in lasso type conveniently.If set to TRUE, this demo aims to use the algorithms in plasso type conveniently(The matrix X must need).Default: FALSE.
#' @param tau Tuning parameter of prior lasso (plasso); using cross-validation can get a great one to calculate. The default is 10.
#' @param meta Take the intersection or union of these methods. (union, intersection) Default: union.
#' @param max.steps Maximum iteration used in calculation. Default: 200.
#'
#' @return A list with the function of parameters ,algorithms, elapse time, grn which a matrix which is the weighted adjacency matrix of the network inferred by the these algorithms.
#' @export
#' @details The demo function aims to call algorithms to infer gene regulatory network.
#' For the existing gene regulatory network algorithms, select several algorithms, and take the intersection or union of the regulatory network results obtained by these methods to obtain a new regulatory network result.
#' @examples
#' A <- matrix(rnorm(200,0,1),20,10)
#' B <- matrix(rnorm(1000,0,1),100)
#' X <- matrix(0,100,20)
#' s <- c(1:10)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' algorithm=c('Pearson','ARACNE.a','CLR','OMP')
#' res <- demo_meta(A=A, B = B,s=s, X=X,algorithm=algorithm,meta='union',max.steps = 200,
#'                  gene=gene,TF = TF,file = NULL,cl.cores=2,verbose = TRUE)
#' grn <- res$grn
demo_meta <- function(A=NULL, B ,s=NULL, X=NULL,algorithm=c('Pearson','ARACNE.a','CLR'),meta='union',max.steps = 200,gene,TF = NULL,p=FALSE,tau=10,file = NULL,cl.cores=2,verbose = TRUE){
  call <- match.call()
  path1 <- getwd()

  if (! is.null(file)) {
    dir <- 'Output'
    if (! dir.exists(dir)) {
      dir.create(dir)
    }
    Outpath <- str_c(dir,'/Output_meta')
    if (! dir.exists(Outpath)) {
      dir.create(Outpath)
    }
    setwd(Outpath)
  }

  colnames(B) <- NULL
  rownames(B) <- gene

  t1 <- Sys.time()
  grn = run_meta(A=A, B = B,s=s, X=X,algorithm=algorithm,meta=meta,max.steps = max.steps,gene = gene,TF = TF, p=p, tau=tau, file = file,cl.cores=cl.cores)
  t2 <- Sys.time()


  if (verbose) {
    print(sprintf("%s has completed", 'meta'))
  }
  time <- as.numeric(t2) - as.numeric(t1)
  if (is.null(file)) {
    object <- list(call = call, algorithm = algorithm, time = time, grn = grn)
  }else{
    object <- list(call = call, algorithm = algorithm, time = time, savePath = Outpath,grn = grn)
  }

  setwd(path1)
  return(object)

}


