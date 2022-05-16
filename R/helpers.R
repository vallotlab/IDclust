#' Create Pseudo-Bulk Matrice from scRNA clusters & replicates
#'
#' @param Seu A Seurat object containing scRNA dataset with 'cell_cluster' 
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
                                      by = "cell_cluster",
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



#' Calculate FDR for SingleCellExperiment object
#'
#' @param scExp A SingleCellExperiment object 
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster.
#' @param limit An integer specifying the minimum number of significantly 
#' enriched / depleted features required in order for a subcluster to be called
#' a 'true' subcluster
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param logFC.th A numeric specifying the log Fold-Change of activation
#' (% cells activated) threshold above which a feature is defined as
#' significantly differential.
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#' @param cluster_of_origin Name of the parent cluster
#' @param nThreads If runFDR==TRUE. An integer specifying of threads to use
#' for the calculation of the FDR.
#' @param iterations An integer specifyung the number of iterations of random
#' subsampling to run.
#' @param ... Additional parameters passed to the differential_function. See 
#' [differential_ChromSCape()] for more information on additional
#' parameters for the default function.
#' 
#' @return
#' @export
#' 
#' @import foreach
#' 
#' @examples
#' if(requireNamespace("ChromSCape", quietly=TRUE)){
#' data("scExp", package = "IDclust")
#' scExp = ChromSCape::find_clusters_louvain_scExp(scExp, resolution = 0.1)
#' FDR_df = calculate_FDR_scEpigenomics(scExp, nThreads = 1, iterations = 2)
#' head(FDR_df)
#' }
calculate_FDR_scEpigenomics <- function(scExp,
                                        differential_function = differential_ChromSCape,
                                        by = "cell_cluster",
                                        limit = 5,
                                        qval.th = 0.01,
                                        logFC.th = log2(2),
                                        min.pct = 0.01,
                                        min_frac_cell_assigned = 0.1,
                                        cluster_of_origin = "Omega",
                                        nThreads = 10,
                                        iterations = 10,
                                        ...){
  if(!requireNamespace("doParallel")){
    stop("IDclust::calculate_FDR_scEpigenomics - Please install ",
         "doParallel to run FDR in parallel")
  } else{
      requireNamespace("doParallel")
  }
  if(!requireNamespace("foreach")){
    stop("IDclust::calculate_FDR_scEpigenomics - Please install ",
         "foreach to run FDR in parallel")
  } else{
      requireNamespace("foreach")
  }
  
  myCluster <- parallel::makeCluster(nThreads, type = "FORK")
  doParallel::registerDoParallel(myCluster)
  FDR_iterations = foreach::foreach(i = 1:iterations) %dopar% {
    set.seed(i * 47)
    scExp$cell_cluster = sample(scExp$cell_cluster, length(scExp$cell_cluster), replace = FALSE)
    DA <- find_differentiated_clusters(object = scExp,
                                       differential_function = differential_function,
                                       by = by,
                                       qval.th = qval.th,
                                       min.pct = min.pct,
                                       logFC.th = logFC.th,
                                       min_frac_cell_assigned = min_frac_cell_assigned,
                                       limit = limit,
                                       limit_by_proportion = NULL,
                                       cluster_of_origin = cluster_of_origin, 
                                       ...)

    df = DA$diffmat_n[,-4]
    df$iteration = i
    df
  }  
  
  FDR = do.call("rbind", FDR_iterations)
  
  parallel::stopCluster(myCluster)
  gc()
  
  return(FDR)
}

#' Retrieve Top differential markers
#'
#' @param scExp  A SingleCellExperiment object 
#' @param IDC_DA_list A list of data.frame of  differential features at each 
#' clustering  iteration produced by [iterative_differential_clustering()]
#' ("IDC_DA.qs").
#' @param top An integer specifying the number of top features to retrieve per
#' cluster.
#' @param gene_col A character specifying the column in which to retrieve the 
#' gene / feature name.
#' @param pseudogene_pattern A character specifying the pattern of 'pseudo-genes'
#' to exclude from the top markers.
#' @param distanceToTSS An integer specifying the maximum distance to consider 
#' a feature close enough to a gene.
#'
#' @return
#' @export
#'
#' @examples
top_differential_markers.SingleCellExperiment <- function(
    scExp, 
    IDC_DA_list,
    top = 1,
    gene_col = "Gene", 
    pseudogene_pattern = "Rik|Vmn|Gm|AW",
    distanceToTSS = 1000
){
  topmarkers = list()
  df = as.data.frame(SingleCellExperiment::rowData(scExp))
  df = df[grep(invert = TRUE, pseudogene_pattern, df[[gene_col]]),]
  
  for(origin in names(IDC_DA_list)){
    top = IDC_DA_list[[origin]]
    top$Gene = df$Gene[match(top$ID, df$ID)]
    top$distanceToTSS = df$distanceToTSS[match(top$ID, df$ID)]
    top$origin = origin
    if(!is.null(distanceToTSS) | !is.na(distanceToTSS))
      top = top[which(top$distanceToTSS < distanceToTSS),] 
    topmarkers[[origin]] = top %>% 
      tidyr::gather(key = "cluster", value = "logFC", starts_with("logFC")) %>% 
      dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = logFC, n =  top, with_ties = FALSE) %>%
      dplyr::select(ID, origin, cluster, logFC, Gene, distanceToTSS)
  }
  df = do.call("rbind",topmarkers)
  df$cluster = paste0(df$origin, ":", gsub("logFC.","",df$cluster))
  return(df)
}