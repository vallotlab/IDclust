#' Normalize and run PCA with Seurat 
#'
#' @description  Runs Seurat PCA with 10 PC and by removing the first 
#' component. See [ChromSCape::preprocessing_filtering_and_reduction()]
#' @param Seu A Seurat object to preprocess with Seurat.
#' @param dim_red The name of the slot to save the dimensionality reduction at
#' each step in the  \code{Seurat::Reductions(Seu)}.
#' @param n_dims  An integer specifying the number of first dimensions to keep 
#' in the dimensionality reduction step.
#'
#' @return A Seurat object with the dimensionality reduction filled
#' @export
#'
#' @examples
#' data(Seu)
#' 
processing_Seurat <- function(Seu, n_dims = 50, dim_red = "pca"){
  Seu = Seurat::FindVariableFeatures(Seu, selection.method = "vst", verbose = FALSE)
  Seu = Seurat::ScaleData(Seu, verbose = FALSE)
  Seu = Seurat::RunPCA(Seu, npcs = n_dims, reduction.name = dim_red,
                       features = Seurat::VariableFeatures(object = Seu), verbose = FALSE)
  return(Seu)
}


#' Normalize and run PCA with ChromSCape (LSI)
#' 
#' @description Runs ChromSCape LSI and  removes the first 
#' component, often driven by total library size.
#'  See [ChromSCape::preprocessing_filtering_and_reduction()].
#' 
#' @param scExp A Seurat object to preprocess with Seurat.
#' @param dim_red The name of the slot to save the dimensionality reduction at
#' each step in the  \code{reducedDimNames(scExp)}.
#' @param n_dims  An integer specifying the number of first dimensions to keep 
#' in the dimensionality reduction step.
#'
#' @return A SingleCellExperiment object with the dimensionality reduction filled
#' @export 
#'
#' @examples
processing_ChromSCape <- function(scExp, n_dims = 10, dim_red = "PCA"){

  # Re-run TFIDF and PCA and clusters using Louvain algorithm 
  scExp = ChromSCape::find_top_features(scExp,
                                        n = nrow(scExp),
                                        keep_others = FALSE)
  scExp = ChromSCape::normalize_scExp(scExp, type = "TFIDF")
  scExp = ChromSCape::reduce_dims_scExp(scExp, dimension_reductions = dim_red,
                                         n = n_dims, remove_PC = "Component_1")
  return(scExp)
}
