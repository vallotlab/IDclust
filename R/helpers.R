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
  raw_mat = Seurat::GetAssayData(object, assay = assay, layer = "counts")
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
          mat[,n] = Matrix::rowSums(raw_mat[,cells, drop = F])
          names_mat = c(names_mat, paste0(i,"_",b))
          rep_done = rep_done + 1
      } else {
          cells_not_used = c(cells_not_used, cells)
      }
      
    }
    if(rep_done < n_rep & rep_done != 0){
      if(length(cells_not_used) > 0){
        n = n + 1
        mat[,n] =  Matrix::rowSums(raw_mat[,cells_not_used, drop = F])
        names_mat = c(names_mat, paste0(i,"_small_rep"))
      }
    }
    if(rep_done == 0){
          n = n + 1
          cells_rep1 = sample(cells_not_used, floor(length(cells_not_used)/2))
          mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep1, drop = F])
          names_mat = c(names_mat, paste0(i,"_small_rep1"))
          
          n = n + 1
          cells_rep2 = cells_not_used[which(! cells_not_used %in% cells_rep1)]
          names_mat = c(names_mat, paste0(i,"_small_rep2"))
          mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep2, drop = F])
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
        mat[,n] = Matrix::rowSums(raw_mat[,cells, drop = F])
        names_mat = c(names_mat, paste0(i,"_",b))
        rep_done = rep_done + 1
      } else {
        cells_not_used = c(cells_not_used, cells)
      }
      
    }
    if(rep_done < n_rep & rep_done != 0){
      n = n + 1
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_not_used, drop = F])
      names_mat = c(names_mat, paste0(i,"_small_rep"))
    }
    if(rep_done == 0){
      n = n + 1
      cells_rep1 = sample(cells_not_used, floor(length(cells_not_used)/2))
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep1, drop = F])
      names_mat = c(names_mat, paste0(i,"_small_rep1"))
      
      n = n + 1
      cells_rep2 = cells_not_used[which(! cells_not_used %in% cells_rep1)]
      names_mat = c(names_mat, paste0(i,"_small_rep2"))
      mat[,n] =  Matrix::rowSums(raw_mat[,cells_rep2, drop = F])
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



#' Reconstruct markers history of IDclusters
#'
#' @description  
#' 
#' First, this function runs a final one vs all differential analysis between
#' all the final IDclusters. The markers obtained this way that are not already
#' present in the IDclust differential analysis are marked as 'global'.
#' Then, for each cluster the markers of the parents that are not already
#' present in the IDclust differential analysis are associated to the given
#' clusters and marked as 'clade-specific'. A new column is added that contains
#' the clade level. The lower the level, the less specific the marker is.
#' Finally, the cluster specific markers are marked as 'cluster-specific' 
#' The rationale is to improve the overall cell type identification and to 
#' return an idea of the specificity of each marker.
#' 
#' 
#' @param object A Seurat or SingleCellObject containing the by metadata
#' column
#' @param IDC_DA The differential analysis data.frame returned by 
#' iterative_differential_clustering function
#' @param by The column containing the clusters ('IDcluster')
#' @param logFC.th The logFC threshold used for global markers.
#' @param qval.th The adjusted p-value thresholds used for global markers
#' @param biological_replicate_col Optional. A character indicating a column
#' containing replicate IDs.
#' @param assay The assay for Seurat object.
#'
#' @return A differential analysis data.frame containing final clusters and for
#' each cluster the specificity of the marker.
#' 
#' @export
#'
#' @examples 
#' data("IDC_DA_scRNA")
#' data("Seu")
#' reconstruct_markers(IDC_DA_scRNA)
reconstruct_markers <- function(object, IDC_DA, by = "IDcluster", 
                                logFC.th = log2(2), qval.th = 0.01,
                                biological_replicate_col = NULL, ...){
  
  res = differential_edgeR_pseudobulk_LRT(object,
                                          by = by,
                                          logFC.th = logFC.th,
                                          qval.th = qval.th,
                                          biological_replicate_col = biological_replicate_col,
                                          ...)
  res$IDcluster = res$cluster
  res$cluster_of_origin = gsub("_[A-Z][0-9]+$","", res$cluster)
  res$cluster = gsub(".*_","", res$cluster)
  res$type = "global"
  res = add_gene_to_DA_list( scExp = object, IDC_DA = res,  feature_ID_col = "ID",
                       gene_col = "Gene",  distanceToTSS = 1000,  split = TRUE,  split_char = ", ")
  IDC_DA$type = "cluster-specific"
  IDC_DA = IDC_DA[order(nchar(IDC_DA$IDcluster), IDC_DA$IDcluster),]
  
  for(cluster_of_origin in unique(IDC_DA$cluster_of_origin)){
    if(cluster_of_origin != "Alpha"){
    clusters = unique(IDC_DA$IDcluster[which(IDC_DA$cluster_of_origin == cluster_of_origin)])
    clade_level_genes = IDC_DA$gene[IDC_DA$IDcluster == cluster_of_origin]
    cluster_level_genes = IDC_DA$gene[IDC_DA$IDcluster %in% clusters]
    clade_level_genes = setdiff(clade_level_genes, cluster_level_genes)
    IDC_DA$gene[IDC_DA$IDcluster == cluster_of_origin & !(IDC_DA$gene %in% clade_level_genes)] = NA
    }
  }
  if(length(which(is.na(IDC_DA$gene)))> 0) IDC_DA = IDC_DA[-which(is.na(IDC_DA$gene)),]
  
  for(cluster in unique(IDC_DA$IDcluster)){
    specific = IDC_DA$gene[IDC_DA$IDcluster == cluster]
    cluster_of_origin = IDC_DA$cluster_of_origin[IDC_DA$IDcluster == cluster][1]
    if(cluster_of_origin != "Alpha"){
      unspecific = IDC_DA[IDC_DA$IDcluster == cluster_of_origin,]
      unspecific = unspecific[!(unspecific$gene %in% specific),]
      unspecific$IDcluster = cluster
      unspecific$cluster_of_origin = cluster_of_origin
      unspecific$type = "clade-specific"
      
      IDC_DA = rbind(IDC_DA, unspecific)
    }
  }
  # IDC_DA = IDC_DA[IDC_DA$IDcluster %in% unlist(object[[by]]),]
  IDC_DA = IDC_DA[,colnames(IDC_DA) %in% colnames(res)]
  res = res[,match(colnames(IDC_DA), colnames(res))]
  unique_features_res = paste0(res$IDcluster,"_", res$gene)
  unique_features_IDC_DA = paste0(IDC_DA$IDcluster,"_", IDC_DA$gene)
  res = res[which(! unique_features_res %in% unique_features_IDC_DA),]
  IDC_DA = rbind(IDC_DA, res)
  IDC_DA = IDC_DA[order(nchar(IDC_DA$IDcluster), IDC_DA$IDcluster),]
  rownames(IDC_DA) = NULL
  return(IDC_DA)
}

#' Find most prevalent cell type in a given list of markers 
#'
#' @param genes A list of genes
#' @param organ_pattern Use to restrict the database with organs matching the 
#' pattern
#' @param species If is not human, use babelgene to translate
#' @param plot_cell_type Plot bar graph for cell types ?
#' @param plot_organs Plot bar graph for organs ?
#'
#' @return
#' @export
#'
#' @examples
#' cell_type_markers(c("ADAMTSL1","ADCY5","ADCY9","AFF3","ALCAM"))
cell_type_markers <- function(genes, organ_pattern = "", species = "human",
                  plot_cell_type = TRUE, plot_organs = FALSE ){
  if(species != "human") {
    try({
    orthologs = babelgene::orthologs(genes, human = FALSE, species = species)
    genes = intersect(orthologs$human_symbol, cell_marker_db$official.gene.symbol)
    } ,silent = T)
  } else{
    genes = intersect(genes, cell_marker_db$official.gene.symbol)  
  }
  
  cell_marker_db. = cell_marker_db[grep(organ_pattern, cell_marker_db$organ, ignore.case = T),]
  total = table(cell_marker_db.$cell.type)
  cell_marker_db. = cell_marker_db.[which(cell_marker_db.$official.gene.symbol %in% genes),]
  
  if(nrow(cell_marker_db.) > 2){

      celltypes = as.data.frame(table(cell_marker_db.$cell.type))
      celltypes$total = total[match(celltypes$Var1, names(total) )]
      
      celltypes = celltypes %>% dplyr::arrange(dplyr::desc(Freq))
      if(plot_cell_type){
      print(celltypes %>% head(20) %>% ggplot(aes(x = forcats::fct_reorder(Var1, Freq), y = Freq)) +
              geom_bar(stat = "identity", fill = "#19297C", position = position_stack(reverse = F)) +
              NoLegend() + themplot + theme(axis.text.x = element_text(angle = 90)) + xlab("") +
              ylab("Number of genes") + geom_text(aes(label=total), position = position_stack(vjust = 1.075)) +
              ggtitle(paste0("Genes = ", length(genes)))
      )
 
      celltypes = celltypes %>% mutate(ratio = Freq / total) %>% dplyr::arrange(dplyr::desc(ratio))

        print(celltypes %>% head(20) %>%  ggplot(aes(x = forcats::fct_reorder(Var1, ratio), y = ratio)) +
                geom_bar(stat = "identity", fill = "#19297C", position = position_stack(reverse = F)) +
                NoLegend() + themplot + theme(axis.text.x = element_text(angle = 90)) + xlab("") +
                ylab("Ratio genes in list / total genes") + geom_text(aes(label=total), position = position_stack(vjust = 1.075)) 
              + ggtitle(paste0("Genes = ", length(genes)))
        )
      }

      cell_marker_db. = cell_marker_db[grep(organ_pattern, cell_marker_db$organ, ignore.case = T),]
      list_celltypes = list()
      for(i in unique(cell_marker_db.$cell.type)){
        list_celltypes[[i]] = cell_marker_db.$official.gene.symbol[cell_marker_db.$cell.type == i]
      }
      
      GeneSetsDf = cell_marker_db.[,c("cell.type","organ")]
      colnames(GeneSetsDf) = c("Gene.Set", "Class")
      enrichments = ChromSCape:::results_enrichmentTest(differentialGenes = genes, 
                                            enrichment_qval = 0.25,
                                            GeneSets = list_celltypes,
                                            GeneSetsDf = GeneSetsDf,
                                            GenePool =  unique(cell_marker_db.$official.gene.symbol))
      enrichments = unique(enrichments)
      enrichments = enrichments %>% group_by(Gene.Set) %>% slice_min(`p-value`, with_ties = F)
      enrichments = enrichments %>% arrange(`p-value`) 
      if(plot_cell_type){
      print(enrichments %>% head(20) %>%  ggplot(aes(x = forcats::fct_reorder(Gene.Set, `p-value`), y = -log10(`p-value`) )) +
              geom_bar(stat = "identity", fill = "#19297C", position = position_stack(reverse = F)) +
              NoLegend() + themplot + theme(axis.text.x = element_text(angle = 90)) + xlab("") +
              ylab("-log10(p-value)") + geom_text(aes(label= Nb_of_deregulated_genes), position = position_stack(vjust = 1.075)) 
            + ggtitle(paste0("Genes = ", length(genes)))
      )
      }
    
    if(plot_organs){
      organs = as.data.frame(table(cell_marker_db.$organ))
      print(
        organs %>% head(20) %>% ggplot(aes(x = forcats::fct_reorder(Var1, Freq), y = Freq)) +
          geom_bar(stat = "identity", fill = "#19297C", position = position_stack(reverse = F)) +
          NoLegend() + themplot + theme(axis.text.x = element_text(angle = 90))
      )
    }
    return(list("cell_types" = celltypes, "enrichment" = enrichments)) 
  }
  
}


