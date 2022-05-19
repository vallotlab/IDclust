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
#' #' # Clustering of scExp scH3K27ac object (Paired-Tag)
#' if(requireNamespace("ChromSCape", quietly=TRUE)){
#'  top_differential_markers.SingleCellExperiment(scExp)
#' data("scExp", package = "IDclust")
#' scExp = iterative_differential_clustering(scExp,  saving = FALSE, plotting =FALSE,
#' logFC.th = 0.5, qval.th = 0.01)
#' 
#' }
#'
#' 
#' 
top_differential_markers.default <- function(
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
top_differential_markers.Seurat <- function(
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