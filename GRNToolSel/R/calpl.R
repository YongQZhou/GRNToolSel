
#' Calculate the predictions of grn and the labels of the high-throughput golden standard.
#'
#' @description Evalute the infer network and high-throughput golden standard matrix. Calculate the predictions of grn and the labels of the high-throughput golden standard.
#' @param HGSmat A matrix which is the high-throughput golden standard. Rows contain genes and columns contain tfs ,the value between in one or zero.
#' @param grn A matrix which is the weighted adjacency matrix of the inferred network by this algorithm or the result from the function of out_grn.
#' @param sym If set to TRUE, only the regulatory relationship in the high-throughput gold standard is considered, and the value of 0 is not considered, so that the number of true positive and false positive is equal, and both false negative and true negative are 0.
#' If set to FALSE, It is assumed that the regulatory relationships that do not exist in the high-throughput gold standard are non regulatory and 0.Default: TRUE.
#' @param num The number of gene regulator imformation. Default: NULL
#'
#' @return A data.frame. The first element, dat$predictions, is a vector of numerical predictions. The second element, dat$labels, is a vector of cordatponding class labels.
#' @export
#'
#' @examples
#' A <- matrix(rnorm(200,0,1),20,10)
#' B <- matrix(rnorm(1000,0,1),100)
#' X <- matrix(0,100,20)
#' s <- c(1:10)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' res <- demo_reg(A, B, s, X, 'ADMM',max.steps = 200, TF, gene, cl.core = 2,file=NULL,verbose=FALSE)
#' grn <- res$grn
#' HGSmat <- matrix(0,100,20)
#' rownames(HGSmat) <- gene
#' colnames(HGSmat) <- TF
#' HGSmat[sample(2000,200)] <- 1
#' dat <- calpl(HGSmat,grn,sym=FALSE)
#' head(dat)
calpl <- function(HGSmat,grn,sym = TRUE,num = NULL){
  if (ncol(grn) <= 3) {
    grnmat = tidyfst::df_mat(grn,grn[,2],grn[,1],grn[,3])
    grnmat <- ifelse(is.na(grnmat),0,grnmat)
    grn <- grnmat
  }
  grn <- as.matrix(grn)
  grn <- abs(grn)
  HGSmat <- abs(HGSmat)
  grn <- grn/max(grn)
  if (!is.null(num)) {
    diag(grn) <- 0
    grn <- tidyfst::mat_df(grn)
    grn2 <- grn[,c(2,1,3)]
    grn2 <- grn2[-which(grn2[,3]==0),]
    if (is.null(num)) {
      num <- nrow(grn2)
    }
    if (num > nrow(grn2)) {
      num <- nrow(grn2)
    }
    grn2 <- grn2[order(abs(grn2[,3]),decreasing = T),]
    grn2 <- grn2[1:num,]
    
    grnmat = tidyfst::df_mat(grn2,grn2[,2],grn2[,1],grn2[,3])
    grnmat <- ifelse(is.na(grnmat),0,grnmat)
    grn <- grnmat
  }
  id1 <- match(rownames(grn),rownames(HGSmat))
  id2 <- match(colnames(grn),colnames(HGSmat))
  HGSmat <- HGSmat[id1,id2]
  
  if (sym) {
    IND1 <- which(HGSmat == 1 & grn != 0, arr.ind = TRUE)
    N <- nrow(IND1)
    if (N != 0) {
      IND <- which(HGSmat == 0 & grn != 0, arr.ind = TRUE)
      N2 <- nrow(IND)
      if (N > N2) {
        N <- N2
      }
      set.seed(123)
      IND <- IND[sample(nrow(IND), N, replace = F), ]
      IND1 <- IND1[sample(nrow(IND1), N, replace = F),
      ]
      IND <- rbind(IND1, IND)
      dat <- grn[IND]
      dat <- matrix(c(dat, c(rep(1, N), rep(0, N))), ncol = 2)
      dat <- as.data.frame(dat)
      colnames(dat) <- c("predictions", "labels")
    }
    else {
      dat <- NULL
    }
  }
  if (!sym) {
    id1 <- which(grn !=0)
    predictions <- grn[id1]
    labels <- HGSmat[id1]
    dat <- data.frame(predictions = as.numeric(predictions),
                      labels = as.numeric(labels))
  }
  return(dat)
}
