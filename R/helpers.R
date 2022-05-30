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
#' mat <- create_pseudobulk_mat_Seu(Seu, by = "seurat_clusters")
#' }
create_pseudobulk_mat_Seu <- function(Seu,
                                      by = "IDcluster",
                                      biological_replicate_col = NULL,
                                      assay = "RNA"){
  cluster_u = unique(Seu@meta.data[[by]])
  if(is.null(biological_replicate_col)){
    Seu$fake_replicate = sample(paste0("rep_",1:3), ncol(Seu), replace = TRUE)
    biological_replicate_col = "fake_replicate"
  }
  biological_replicates = unique(unlist(Seu[[biological_replicate_col]]))
  
  mat = matrix(0, nrow = nrow(Seu@assays[[assay]]@counts), ncol = length(cluster_u) * length(biological_replicates))
  rownames(mat) = rownames(Seu@assays[[assay]]@counts)
  
  n = 0
  names_mat =c()
  for(i in cluster_u){
    for(b in biological_replicates){
      n = n+1
      cells = colnames(Seu)[which(Seu@meta.data[[by]] == i & unlist(Seu[[biological_replicate_col]]) == b)]
      if(length(cells) > 25) {
        mat[,n] = Matrix::rowSums(Seu@assays[[assay]]@counts[,cells])
        names_mat = c(names_mat,paste0(i,"_",b))
      }
    }
  }
  mat = mat[,which(Matrix::colSums(mat) > 0)]
  colnames(mat) = names_mat
  return(mat)
}



