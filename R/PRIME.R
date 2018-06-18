
#' Probabilistic imputation to reduce dropout effects in single cell sequencing data
#'
#' Probabilistic imputation to reduce dropout effects in single cell sequencing data
#'
#'
#' @param sc_cnt (M x N) dimensional input matrix. This matrix sould be in count scale (not a log scale)
#'  M: number of genes and N: Number of cells
#' @param max_it Maximum number of iteration
#' @param alp Shape parameter for a sigmoid function (probabilistic weight)
#' @param err_max Error tolerance to stop the iteration
#' @return Normalized imputation result in count scale
#'
#' @author Hyundoo Jeong
#'
#' @references
#'  Hyundoo Jeong and Zhandong Liu
#'
#'  PRIME: a probabilistic imputation method to reduce dropout effects in single cell RNA sequencing
#'
#' @examples
#' data(testdata)
#' PRIME_res <- PRIME(testdata)
#' pca_res <-prcomp(t(log10(1+PRIME_res)))
#' plot(pca_res$x[,1], pca_res$x[,2], bg = c("red", "blue","green")[factor(label)],  type = "p", pch = 21, xlab = "PC1", ylab = "PC2")
#' legend("topleft", title="Cell types", c("G1","S","G2M"), fill = c("red", "blue","green"), horiz = TRUE)
#'
PRIME <- function(sc_cnt, max_it = 5, alp = 1, err_max = 0.05){
  message('Start PRIME')
  set.seed(1234)
  # normalization
  sc_data_raw <- log10(1+cpm(sc_cnt))

  sc_data_imputed <- sc_data_raw
  nmse <- Inf
  iter <- 0
  while((nmse > err_max) & (iter < max_it)){
    iter <- iter + 1
    message(paste(iter, '-th iteration ', sep = ''))
    sc_data_prev <- sc_data_imputed

    seurat <- CreateSeuratObject(raw.data = sc_data_prev, min.cells = 0, min.genes = 0)
    seurat <- FindVariableGenes(object = seurat, mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 1, x.high.cutoff = 10, y.cutoff = 1, do.plot = F, display.progress = F)
    var_genes <- seurat@var.genes

    sel_data <- sc_data_prev[var_genes,]
    sel_data <- apply(sel_data, 2, as.numeric)

    my_pca_log <- prcomp_irlba(t(sc_data_prev[var_genes,]), n = 20)
    PCs <- t(my_pca_log$x)
    cor_dist <- cor(PCs, method ='pearson')
    cdist_mat <- FindNBR(cor_dist, qth = 0.9, th = 0.85)

    edist <- EucDist(PCs)
    edist_mat <- FindNBR(as.matrix(edist), qth = 0.9, th = 0.85)

    sim_mat <- 0.5*(edist_mat + cdist_mat)
    sim_mat <- sim_mat^2
    D1 <- Matrix(diag(1/colSums(sim_mat)), sparse = TRUE)
    E1 <- sim_mat %*% D1
    E2 <- (E1 %*% E1)

    sc_data_imputed <- prime_core(as.matrix(E2), as.matrix(sim_mat), sc_data_prev, alp)
    rownames(sc_data_imputed) <- rownames(sc_data_raw)
    colnames(sc_data_imputed) <- colnames(sc_data_raw)

    nmse <- sum(((sc_data_prev - sc_data_imputed)^2))/prod(dim(sc_data_imputed))
  }


  sc_cnt_imputed <- 10^sc_data_imputed - 1
  sc_data_final <- cpm(sc_cnt_imputed)

  message('Done!')
  return(sc_data_final)
}







#' Compute Euclidean similarity
#'
#' Compute Euclidean similarity using a Gaussian kernel.
#' First, it normalize the input
#' then it transforms the distance to similarity using a Gaussian kernel
#' so that similar objects have larger values
#'
#'
#' @param data (M x N) dimensional inputmatrix. N is the number of objects M is the number of features describing each object.
#' @return (N x N) dimensional normalized Euclidean similarity matrix.
#'
#' @author Hyundoo Jeong
#'
#' @references
#'  Hyundoo Jeong and Zhandong Liu
#'
#'  PRIME: a probabilistic imputation method to reduce dropout effects in single cell RNA sequencing
#'
#' @examples
#' EucDist(matrix)
#'
EucDist <- function(data){
  e_dist <- parDist(t(data), method = 'euclidean')
  e_dist <- as.matrix(e_dist)

  max_dist <- max(e_dist)
  t <- 2*(e_dist/max_dist)
  edist.t <- exp(-t^2)

  return(edist.t)
}



#' Find neighboring nodes in the network
#'
#' Find neighboring nodes in the network
#'
#' @param dist_mat Distance matrix
#' @param qth Threshold for the quantile of the similarity
#' @param th Hard threshold to select the neighborind nodes
#' @return Adjacency matrix for the correspondence network
#'
#' @author Hyundoo Jeong
#'
#' @references
#'  Hyundoo Jeong and Zhandong Liu
#'
#'  PRIME: a probabilistic imputation method to reduce dropout effects in single cell RNA sequencing
#'
#' @examples
#' FindNBR(matrix, qth = 0.9, th = 0.85)
#'
FindNBR <- function(dist_mat, qth = 0.9, th = 0.85){
  row_q <- rowQuantiles(dist_mat, probs = qth)
  my_th <- apply(as.matrix(row_q), 1, min, th)
  submat <- apply(dist_mat, 2,  function(x) x - my_th > 0)
  dir_net <- apply(submat, 2, as.numeric)

  adj_mat <- 0.5*(dir_net + t(dir_net))
  adj_mat <- Matrix(adj_mat, sparse = TRUE)
  return(adj_mat)
}
