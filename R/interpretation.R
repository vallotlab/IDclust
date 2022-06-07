
#'  Add genes to differential analysis using a SingleCellExperiment 
#'  
#' @param scExp  A SingleCellExperiment object. The rowData must contain a 
#' column with genes (see param 'gene_col').  
#' @param IDC_DA A data.frame of differential features at each 
#' clustering  iteration produced by [iterative_differential_clustering()]
#' ("IDC_DA.csv").
#' @param feature_ID_col A character specifying the column in which to retrieve 
#' the feature ID.
#' @param gene_col A character specifying the column in which to retrieve the 
#' gene / feature name.
#' @param distanceToTSS An integer specifying the maximum distance to consider 
#' a feature close enough to a gene.
#' @param split A logical indicating wether to split the gene column or not
#' @param split_char If split is TRUE. A character by which the gene column is 
#' splitted.
#' 
#' @return A data.frame with the gene column filled. Might contain more rows
#' than the original data.frame if multiple genes per feature were contained in
#' the rowData of the SingleCellExperiment. 
#' 
#' @export
#' @rdname add_gene_to_DA_list
#' @examples 
#' data("scExp", package = "IDclust")
#' data("IDC_DA_scEpigenomics", package = "IDclust")
#' 
#' IDC_DA_scEpigenomics = add_gene_to_DA_list(
#'     scExp = scExp, 
#'     IDC_DA = IDC_DA_scEpigenomics,
#'     feature_ID_col = "ID", 
#'     gene_col = "Gene", 
#'     distanceToTSS = 1000,
#'     split = TRUE,
#'     split_char = ", "
#' )
#' 
add_gene_to_DA_list <- function(
    scExp, 
    IDC_DA,
    feature_ID_col = "ID", 
    gene_col = "Gene", 
    distanceToTSS = 1000,
    split = TRUE,
    split_char = ", "
){
    df = as.data.frame(SingleCellExperiment::rowData(scExp))

    IDC_DA_list = list()
    for(origin in unique(IDC_DA$cluster_of_origin)){
        topmarker = IDC_DA %>% dplyr::filter(cluster_of_origin == origin)
        topmarker[[gene_col]] = df[[gene_col]][match(topmarker[[feature_ID_col]], df[[feature_ID_col]])]
        if(!is.null(distanceToTSS)) topmarker$distanceToTSS = df$distanceToTSS[match(topmarker[[feature_ID_col]], df[[feature_ID_col]])]
        if(!is.null(distanceToTSS) | !is.na(distanceToTSS))
            topmarker = topmarker[which(topmarker$distanceToTSS < distanceToTSS),] 
        if(split){
          topmarker = topmarker %>% tidyr::separate_rows(.data[[gene_col]], sep = split_char )
        }
        IDC_DA_list[[origin]] = topmarker
    }
    IDC_DA = do.call("rbind", IDC_DA_list)
    rownames(IDC_DA) = NULL
    
    return(IDC_DA)
}




#' Retrieve top marker genes from list of differential analysis
#'
#' @param object  A data.frame of  differential features at each 
#' clustering  iteration produced by [iterative_differential_clustering()]
#' ("IDC_DA.csv").
#' @param top An integer specifying the number of top features to retrieve per
#' cluster.
#' @param gene_col A character specifying the column in which to retrieve the 
#' gene / feature name.
#' @param logFC_col A character specifying the column in which to retrieve the 
#' logFC.
#' @param qvalue_col A character specifying the column in which to retrieve the 
#' adjusted p.value.
#' @param order_by A character specifying the column by which to order the top
#' markers (default to logFC_col.
#' @param pseudogene_pattern A character specifying the pattern of 'pseudo-genes'
#' to exclude from the top markers.
#' 
#' @return
#' @export
#' @rdname top_differential_markers
#' @examples 
#' #scRNA
#' data("IDC_DA_scRNA", package = "IDclust")
#'
#' top_differential_markers(
#'     IDC_DA_scRNA,
#'     top = 1,
#'     gene_col = "gene",
#'     logFC_col = "avg_log2FC",
#'     qvalue_col = "p_val_adj",
#'     order_by = "logFC_col",
#'     pseudogene_pattern = NULL
#' )
#'
#' #scEpigenomics
#' data("scExp", package = "IDclust")
#' data("IDC_DA_scEpigenomics", package = "IDclust")
#' 
#' # We must first add the gene information to the DA list: 
#' IDC_DA_scEpigenomics = add_gene_to_DA_list(
#'     scExp = scExp, 
#'     IDC_DA = IDC_DA_scEpigenomics,
#'     feature_ID_col = "ID", 
#'     gene_col = "Gene", 
#'     distanceToTSS = 1000,
#'     split = TRUE,
#'     split_char = ", "
#' )
#' 
#' top_differential_markers(
#'     IDC_DA_scEpigenomics,
#'     top = 3,
#'     gene_col = "Gene",
#'     logFC_col = "logFC",
#'     qvalue_col = "qval",
#'     order_by = "logFC_col",
#'     pseudogene_pattern = "Rik|Vmn|Gm|AW"
#' )
#'     
top_differential_markers <- function(
    object,
    top = 1,
    gene_col = "gene",
    logFC_col = "avg_log2FC",
    qvalue_col = "p_val_adj",
    order_by = c("logFC_col", "qvalue_col")[1],
    pseudogene_pattern = NULL
){
    topmarkers = list()
    if(order_by == "logFC_col") order_by = logFC_col else order_by = qvalue_col
    for(origin in unique(object$cluster_of_origin)){
        topmarker = object %>% dplyr::filter(cluster_of_origin == origin)

        if(!is.null(pseudogene_pattern)) topmarker = topmarker[grep(pseudogene_pattern, topmarker[[gene_col]], invert = TRUE),]
        
        topmarkers[[origin]] = topmarker %>%
            dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = .data[[order_by]], n =  top, with_ties = FALSE) %>%
            dplyr::select(cluster_of_origin, cluster, .data[[logFC_col]], .data[[qvalue_col]], .data[[gene_col]])
    }
    
    df = do.call("rbind",topmarkers)
    rownames(topmarkers) = NULL
    
    return(df)
}

#'  Find enriched pathways in marker genes of IDclust 
#'  
#' @param IDC_DA A data.frame of  differential features at each 
#' clustering  iteration produced by [iterative_differential_clustering()]
#' ("IDC_DA.qs").
#' @param top An integer specifying the number of top pathways to retrieve per
#' cluster.
#' @param gene_col A character specifying the column in which to retrieve the 
#' gene / feature name.
#' @param qval.th A numeric specifying the adjusted p.value below which a 
#' pathway is considered as significantly enriched
#' @param website site requested
#' @param databases (Required). Character vector of databases to search.
#'  See https://maayanlab.cloud/Enrichr/ for available databases.
#' @param order_by_database A logical. If TRUE, row will appear in the order of
#' the databases vector, then within each database are sorted by adjusted
#' pvalue.
#' @param max_genes_enriched An integer specifying the maximum number of genes
#' per cluster to enrich. If there are more than max_genes_enriched genes in the
#' list, will keep only the top max_genes_enriched genes for enrichment.
#' 
#' @return
#' @export
#' @rdname top_enriched_pathways
#' 
#' @examples
#' if(requireNamespace("enrichR")){
#' #scRNA
#' data("IDC_DA_scRNA", package = "IDclust")
#' library(enrichR)
#' top_enriched_pathways(
#'     IDC_DA_scRNA,
#'     top = 1,
#'     gene_col = "gene",
#'     qval.th = 0.1)
#'
#' #scEpigenomics
#' data("scExp", package = "IDclust")
#' data("IDC_DA_scEpigenomics", package = "IDclust")
#' 
#' # We must first add the gene information to the DA list: 
#' IDC_DA_scEpigenomics = add_gene_to_DA_list(
#'     scExp = scExp, 
#'     IDC_DA = IDC_DA_scEpigenomics,
#'     feature_ID_col = "ID", 
#'     gene_col = "Gene", 
#'     distanceToTSS = 1000,
#'     split = TRUE,
#'     split_char = ", "
#' )
#' 
#' top_enriched_pathways(
#'     IDC_DA_scEpigenomics,
#'     top = 1,
#'     gene_col = "Gene",
#'     qval.th = 0.1
#' )
#' }
top_enriched_pathways <- function(
    IDC_DA,
    top = 1,
    qval.th = 0.1,
    gene_col = "gene",
    website = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")[1],
    databases = c("MSigDB_Hallmark_2020", "GO_Molecular_Function_2021", "MSigDB_Oncogenic_Signatures"),
    order_by_database = TRUE,
    max_genes_enriched = 1000
){
  enriched_list = list()
  
  for(origin in unique(IDC_DA$cluster_of_origin)){
    topmarker_all = IDC_DA %>% dplyr::filter(cluster_of_origin == origin)
    
    for(clust in unique(topmarker_all$cluster)){
      name = paste0(origin, ":", clust)
      
      topmarker = topmarker_all %>% dplyr::filter(cluster == clust)
      cat("Enriching", nrow(topmarker), "marker genes from cluster", name,"...\n")
      
      genes = topmarker[[gene_col]]
      if(length(genes) > max_genes_enriched){
          
          cat("More than", max_genes_enriched, " genes in the list, only ",
              "keeping the first", max_genes_enriched, "genes in the list.\n")
          genes = head(genes, n = max_genes_enriched)
      } 
          
      enriched_path = differential_pathway(
        genes = genes, 
        qval.th = qval.th, 
        website = website,
        databases = databases,
        order_by_database = order_by_database
      )
      
      if(nrow(enriched_path) < 1){
        
      } else{
          enriched_path$cluster = clust
          enriched_path$cluster_of_origin = origin
        enriched_list[[name]] = enriched_path[1:min(top, nrow(enriched_path)),]
      }
    }
  }
  
  enriched_df = do.call("rbind",enriched_list)
  rownames(enriched_df) = NULL
  
  return(enriched_df)
}

#'  Find differential pathways with enrichR
#'  
#' @param genes A character vector containing at least 10 genes
#' @param qval.th A numeric specifying the adjusted p.value below which a 
#' pathway is considered as significantly enriched
#' @param website site requested
#' @param databases (Required). Character vector of databases to search.
#'  See https://maayanlab.cloud/Enrichr/ for available databases.
#' @param order_by_database A logical. If TRUE, row will appear in the order of
#' the databases vector, then within each database are sorted by adjusted
#' pvalue.
#' 
#' @return A dataframe with enriched pathways.
#' @export
#' @rdname differential_pathway
#' 
#' @import enrichR
#' 
#' @examples
#' 
#' 
#' pathways <- differential_pathway(
#' c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"))
#' 
#' @references Chen, E.Y., Tan, C.M., Kou, Y. et al. Enrichr: interactive and 
#' collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 
#' 128 (2013). https://doi.org/10.1186/1471-2105-14-128
differential_pathway <- function(
    genes,
    qval.th = 0.1,
    website = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")[1],
    databases = c("MSigDB_Hallmark_2020", "GO_Molecular_Function_2021",  "MSigDB_Oncogenic_Signatures"),
    order_by_database = TRUE
){
    
    if(requireNamespace("enrichR")){
    if(length(genes) >= 5){
        enrichR::setEnrichrSite(website) # Human genes
        enriched <- enrichR::enrichr(genes, databases)
        for(i in seq_along(enriched)){
          if(nrow(enriched[[i]] > 0)) enriched[[i]]$database = names(enriched)[i]
        }
        enriched = do.call("rbind", enriched)
        enriched = enriched %>% dplyr::arrange(Adjusted.P.value)
        enriched = enriched %>% dplyr::filter(Adjusted.P.value < qval.th)
        if(order_by_database) enriched = enriched[order(match(enriched$database,databases)),]
        
        return(enriched)
    } else{
        message("Not enough genes (< 5) in the gene list, skipping...")
      return(data.frame())
    }
    } else{
        message("For this function, you need the enrichR package. Install it with install.packages('enrichR')")
        return(data.frame())
    }
}
