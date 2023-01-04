#' Create Pseudo-Bulk Matrice from scRNA clusters & replicates
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @rdname create_pseudobulk_mat
#' @export create_pseudobulk_mat
#' 
#' @examples
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' data("Seu", package = "IDclust")
#' mat <- create_pseudobulk_mat(Seu, by = "seurat_clusters")
#' }
#' 
#' #' if(requireNamespace("SingleCellExperiment", quietly=TRUE)){
#' data("scExp", package = "IDclust")
#' mat <- create_pseudobulk_mat(scExp, by = "cell_clusters")
#' }
create_pseudobulk_mat <- function(object, ...) {
  UseMethod(generic = 'create_pseudobulk_mat', object = object)
  
}

#' Create Pseudo-Bulk Matrice from scRNA clusters & replicates
#'
#' @param Seu A Seurat object containing scRNA dataset with 'IDcluster' 
#' column.
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param biological_replicate_col Optional. A column of the Seurat object 
#' indicating the replicates or batches of the dataset in order to take in 
#' account biological/technical noise. If NULL, will create random layers of 
#' fake replicates.
#' @param assay Assay to use.
#' 
#' @return A pseudo-bulk matrice of cluster spread by replicates / batches /
#' fake replicates.
#' 
#' @export
#' @importFrom Matrix colSums rowSums
#' 
#' @examples
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' data("Seu", package = "IDclust")
#' mat <- create_pseudobulk_mat(Seu, by = "seurat_clusters")
#' }
#' 
create_pseudobulk_mat.Seurat <- function(object,
                                      by = "IDcluster",
                                      biological_replicate_col = NULL,
                                      assay = "RNA"){
  raw_mat = object@assays[[assay]]@counts
  meta = object@meta.data
  
  cluster_u = unique(meta[[by]])
  if(is.null(biological_replicate_col)){
    object$fake_replicate = sample(paste0("rep_",1:3), ncol(object), replace = TRUE)
    biological_replicate_col = "fake_replicate"
  }
  biological_replicates = unique(unlist(object[[biological_replicate_col]]))
  n_rep = length(biological_replicates)
  mat = matrix(0, nrow = nrow(raw_mat), ncol = length(cluster_u) * length(biological_replicates))
  rownames(mat) = rownames(raw_mat)
  
  n = 0
  names_mat =c()
  for(i in cluster_u){
      rep_done = 0
      cells_not_used = c()
    for(b in biological_replicates){
      cells = colnames(object)[which(meta[[by]] == i & unlist(object[[biological_replicate_col]]) == b)]
      if(length(cells) > 25) {
          n = n+1
          mat[,n] = Matrix::rowSums(raw_mat[,cells])
          names_mat = c(names_mat, paste0(i,"_",b))
          rep_done = rep_done + 1
      } else {
          cells_not_used = c(cells_not_used, cells)
      }
      
    }
    if(rep_done < n_rep & rep_done != 0){
        n = n + 1
        mat[,n] =  Matrix::rowSums(raw_mat[,cells_not_used])
        names_mat = c(names_mat, paste0(i,"_small_rep"))
    }
    if(rep_done == 0){
          n = n + 1
          cells_rep1 = sample(cells_not_used, floor(length(cells_not_used)/2))
          mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep1])
          names_mat = c(names_mat, paste0(i,"_small_rep1"))
          
          n = n + 1
          cells_rep2 = cells_not_used[which(! cells_not_used %in% cells_rep1)]
          names_mat = c(names_mat, paste0(i,"_small_rep2"))
          mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep2])
      }
  }
  mat = mat[,which(Matrix::colSums(mat) > 0)]
  colnames(mat) = names_mat
  return(mat)
}

#' Create Pseudo-Bulk Matrix from scEpigenomics clusters & replicates
#'
#' @param object A SingleCellExperiment object containing scEpigenomics dataset 
#' with 'IDcluster' column.
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param biological_replicate_col Optional. A column of the 
#' SingleCellExperiment object indicating the replicates or batches of the
#'  dataset in order to take in account biological/technical noise. If NULL,
#'  will create random layers of fake replicates.
#' 
#' @return A pseudo-bulk matrix of cluster spread by replicates / batches /
#' fake replicates.
#' 
#' @export
#' @importFrom Matrix colSums rowSums
#' 
#' @examples
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' data("Seu", package = "IDclust")
#' mat <- create_pseudobulk_mat.default(object, by = "seurat_clusters")
#' }
create_pseudobulk_mat.default <- function(object,
                                          by = "IDcluster",
                                          biological_replicate_col = NULL){
  
  raw_mat = SingleCellExperiment::counts(object)
  meta = SingleCellExperiment::colData(object)
  
  cluster_u = unique(meta[[by]])
  
  if(is.null(biological_replicate_col)){
    object$fake_replicate = sample(paste0("rep_",1:3), ncol(object), replace = TRUE)
    biological_replicate_col = "fake_replicate"
  }
  biological_replicates = unique(unlist(object[[biological_replicate_col]]))
  n_rep = length(biological_replicates)
  
  mat = matrix(0, nrow = nrow(raw_mat),
               ncol = length(cluster_u) * length(biological_replicates))
  rownames(mat) = rownames(raw_mat)
  
  n = 0
  names_mat =c()
  for(i in cluster_u){
    rep_done = 0
    cells_not_used = c()
    for(b in biological_replicates){
      cells = colnames(object)[which(meta[[by]] == i & unlist(object[[biological_replicate_col]]) == b)]
      if(length(cells) > 25) {
        n = n+1
        mat[,n] = Matrix::rowSums(raw_mat[,cells])
        names_mat = c(names_mat, paste0(i,"_",b))
        rep_done = rep_done + 1
      } else {
        cells_not_used = c(cells_not_used, cells)
      }
      
    }
    if(rep_done < n_rep & rep_done != 0){
      n = n + 1
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_not_used])
      names_mat = c(names_mat, paste0(i,"_small_rep"))
    }
    if(rep_done == 0){
      n = n + 1
      cells_rep1 = sample(cells_not_used, floor(length(cells_not_used)/2))
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep1])
      names_mat = c(names_mat, paste0(i,"_small_rep1"))
      
      n = n + 1
      cells_rep2 = cells_not_used[which(! cells_not_used %in% cells_rep1)]
      names_mat = c(names_mat, paste0(i,"_small_rep2"))
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep2])
    }
  }
  mat = mat[,which(Matrix::colSums(mat) > 0)]
  colnames(mat) = names_mat
  return(mat)
}

#' Summarise Differential Analysis table
#'
#' @param res A differential analysis data.frame
#' @param chr_col Character specifying the chromosome column.
#' @param start_col Character specifying the start column.
#' @param end_col Character specifying the end column.
#' @param ID_col Character specifying the ID column.
#'
#' @return A gathered data.frame with a cluster column
#' @export
#'
#' @examples
summarise_DA <- function(res,
                         chr_col = "chr",
                         start_col = "start",
                         end_col = "end",
                         ID_col = "ID"
){
    res = res %>% dplyr::select(-.data[[chr_col]], -.data[[start_col]], -.data[[end_col]])
    res = res %>% tidyr::gather("key", "var", -.data[[ID_col]])
    res = res %>% tidyr::extract(key, c("column", "cluster"), "(.*)\\.(.*)")
    res = res %>% tidyr::spread(column, var)
    return(res)
}


