#' Output grn
#' @description Transfer matrix data converted to three columns
#' @param grn A matrix which is the weighted adjacency matrix of the inferred network by this algorithm.
#' @param num The number of gene regulator imformation. Default: NULL.
#' @param file A connection, or a character string naming the file to write to. The default not save. Default: NULL.
#' @param TF The Transcription Factors (TFs) The default value NULL means that all the genes are used as candidate regulators. Default: NULL.
#' @param gene A vector of characters containing the target genes.
#' @details The output is consistent with the data format required by Cytoscape. This output data can be used to map the gene regulatory network with Cytoscape.

#' @return A data.frame. The first column is TF, the second is gene, and the last columns is the weighted adjacency values.
#' @export
#'
#' @examples
#' B <- matrix(rnorm(1000,0,1),100)
#' TF <- paste0('G',1:20)
#' gene <- paste0('G',1:100)
#' grn <- run_aracne.a(B,gene,TF)
#' net <- out_grn(grn,TF=TF,gene=gene)
out_grn <- function(grn,num = NULL,TF,gene,file=NULL){
  grn <- as.matrix(grn)
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
  colnames(grn2) <- c('TF','Gene','Weight')
  return(grn2)
}
