#' Differential function Seurat (wilcoxon)
#'
#' @description  Run [Seurat::FindAllMarkers()]  function on the
#' clusters in the 'by' column and return a data.frame  containing the marker 
#' genes of each cluster passing the thresholds .
#' 
#' @param Seu A Seurat object containing scRNA dataset with a metadata column
#'  name matching the by parameter
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#' @param ... Additional parameters to pass to [Seurat::FindAllMarkers()].
#'
#' @seealso See [Seurat::FindAllMarkers()]
#' 
#' @return
#' @export
#'
#' @examples
differential_Seurat <- function(Seu,
                                by = "cell_cluster",
                                logFC.th = log2(1.5),
                                qval.th = 0.01,
                                min.pct = 0.1,
                                ...){
  Seurat::Idents(Seu) = unique(Seu@meta.data[,by])
  res =  Seurat::FindAllMarkers(Seu,
                                logfc.threshold = logFC.th,
                                return.thresh = qval.th,
                                min.pct = min.pct, 
                                only.pos = TRUE,
                                ...)
  return(res)
}

#' Pseudo-bulk differential analysis with edgeR (LRT) 
#'
#' @details Concatenate single-cells into replicates by cluster in order to 
#' create a 'pseudo-bulk' matrice of multiple replicates per cluster. If no 
#' replicates are present, will assign replicates at random to create 3 replicates
#' per cluster. Conducts 'LRT' (likelihood ratio tests) edgeR tests to test. 
#' See [edgeR::glmLRT()] 
#' 
#' @param Seu A Seurat object containing scRNA dataset with a metadata column
#'  name matching the by parameter
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param assay Assay to use.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#' @seealso See [edgeR::glmLRT()]
#'
#' @return A data.frame containing the results of the differential analysis 
#' 
#' @export
#' 
#' @import dplyr
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmFit
#' glmLRT
#' @importFrom Matrix Matrix rowSums
#' @importFrom stats p.adjust
#' 
#' @examples
differential_edgeR_pseudobulk_LRT <- function(Seu,
                                              by = "cell_cluster",
                                              assay = "RNA",
                                              biological_replicate_col = NULL,
                                              logFC.th = log2(1.5),
                                              qval.th = 0.01,
                                              min.pct = 0.1
){
  cluster_u = unique(Seu@meta.data[,by])
  n_cell_assigned = 0
  mat = create_pseudobulk_mat_Seu(Seu, biological_replicate_col, assay)
  res = data.frame("p_val" = 0, "avg_log2FC"= 0, "pct.1" = 0, "pct.2"= 0, "p_val_adj" = 0, "cluster"= "", "gene"= "")
  
  for(i in seq_along(cluster_u)){
    
    group = rep(1, ncol(mat))
    group[grep(paste0("^",cluster_u[i],"_"), colnames(mat))] = 2
    
    if(length(grep(paste0("^",cluster_u[i],"_"), colnames(mat))) > 1 &
       length(grep(paste0("^",cluster_u[i],"_"), colnames(mat), invert = T)) > 1){
      
      group <- as.factor(group)
      y <- edgeR::DGEList(counts=mat, group=group)
      keep <- edgeR::filterByExpr(y)
      if(length(which(keep)) > 0){
        y <- y[keep,,keep.lib.sizes=FALSE]
        y <- edgeR::calcNormFactors(y)
        design <- model.matrix(~group)
        y <- edgeR::estimateDisp(y,design)
        fit <- edgeR::glmFit(y,design)
        lrt <- edgeR::glmLRT(fit,coef=2)
        tab = lrt$table
        
        binmat = Matrix::Matrix((Seu@assays[[assay]]@counts > 0) + 0, sparse = TRUE)
        pct.1 = Matrix::rowSums(binmat[,which(Seu@meta.data[,by] == cluster_u[i])]) / length(which(Seu@meta.data[,by] == cluster_u[i]))
        pct.2 = Matrix::rowSums(binmat[,which(Seu@meta.data[,by]!= cluster_u[i])]) / length(which(Seu@meta.data[,by] != cluster_u[i]))
        
        tab = tab %>% dplyr::filter(abs(logFC) > 0.1 & PValue < 0.1) # very loose filter
        
        res. = data.frame(
          "p_val" = tab$PValue,
          "avg_log2FC"= tab$logFC,
          "pct.1" = pct.1[match(rownames(tab), names(pct.1))], 
          "pct.2"= pct.2[match(rownames(tab), names(pct.2))],
          "p_val_adj" = stats::p.adjust(tab$PValue, method = "bonferroni"),
          "cluster"= cluster_u[i],
          "gene"= rownames(tab)
        )
        res = rbind(res, res.)
      }
    } else{
      cat("Not enough cells to form 2 replicates ... assigning 0 differential genes.\n")
    }
  }
  res = res[-1,]
  
  res = res[which(res$avg_log2FC > logFC.th &
                    res$p_val_adj < qval.th &
                    res$pct.1 > min.pct),]
  
  return(res)
}


#' Differential analysis with ChromSCape (differential activation)
#'
#' Runs [ChromSCape::differential_activation()] function on clusters in the 
#' 'by' column of the SingleCellExperiment object and returns a data.frame 
#' containing the differential features that passed the thresholds.
#' 
#' @param scExp A SingleCellExperiment object containing scRNA dataset with a
#'  metadata column name matching the by parameter
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#' 
#' @seealso See [ChromSCape::differential_activation()]
#' @return
#' @export
#'
#' @examples
differential_ChromSCape <- function(
    scExp,
    by = "cell_cluster",
    logFC.th = log2(1.5),
    qval.th = 0.01,
    min.pct = 0.01
    ){
  
  res = ChromSCape::differential_activation(scExp = scExp, group_by = by)
  res = res %>% dplyr::select(-chr, -start, -end)
  res = res %>% tidyr::gather("key", "var", -ID)
  res = res %>% tidyr::extract(key, c("column", "cluster"), "(.*)\\.(.*)")
  res = res %>% tidyr::spread(column, var)

  res = res %>% dplyr::filter(logFC > logFC.th  & 
                                qval < qval.th  & 
                                group_activation > min.pct)
  
  return(res)
}



