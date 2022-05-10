#' Pseudo-bulk differential analysis (LRT) edgeR 
#'
#' @details Concatenate single-cells into replicates by cluster in order to 
#' create a 'pseudo-bulk' matrice of multiple replicates per cluster. If no 
#' replicates are present, will assign replicates at random to create 3 replicates
#' per cluster. Conducts 'LRT' (likelihood ratio tests) edgeR tests to test 
#' 
#' @param Seu A Seurat object containing scRNA dataset with 'seurat_clusters' 
#' column
#' @param biological_replicate_col Optional. A column of the Seurat object 
#' indicating the replicates or batches of the dataset in order to take in 
#' account biological/technical noise. If NULL, will create random layers of 
#' fake replicates.
#' @seealso glmFit {edgeR}
#'
#' @return A data.frame containing the results of the differential analysis 
#' 
#' @export
#' 
#' @import dplyr
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmFit
#' @importFrom Matrix Matrix rowSums
#' @importFrom stats p.adjust
#' 
#' @examples
pseudobulk_LRT_edgeR <- function(Seu, biological_replicate_col = NULL){
    cluster_u = unique(Seu$seurat_clusters)
    n_cell_assigned = 0
    mat = create_pseudobulk_mat_Seu(Seu, biological_replicate_col)
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
                
                binmat = Matrix::Matrix((Seu@assays$RNA@counts > 0) + 0, sparse = TRUE)
                pct.1 = Matrix::rowSums(binmat[,which(Seu$seurat_clusters == cluster_u[i])]) / length(which(Seu$seurat_clusters == cluster_u[i]))
                pct.2 = Matrix::rowSums(binmat[,which(Seu$seurat_clusters != cluster_u[i])]) / length(which(Seu$seurat_clusters != cluster_u[i]))
                
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
    return(res)
}

#' Create Pseudo-Bulk Matrice from scRNA clusters & replicates
#'
#' @param Seu A Seurat object containing scRNA dataset with 'seurat_clusters' 
#' column.
#' @param biological_replicate_col Optional. A column of the Seurat object 
#' indicating the replicates or batches of the dataset in order to take in 
#' account biological/technical noise. If NULL, will create random layers of 
#' fake replicates.
#' 
#' @return A pseudo-bulk matrice of cluster spread by replicates / batches /
#' fake replicates.
#' 
#' @export
#' @importFrom Matrix colSums rowSums
#' 
#' @examples
create_pseudobulk_mat_Seu <- function(Seu, biological_replicate_col){
    cluster_u = unique(Seu$seurat_clusters)
    if(is.null(biological_replicate_col)){
        Seu$fake_replicate = sample(paste0("rep_",1:3), ncol(Seu), replace = T)
        biological_replicate_col = "fake_replicate"
    }
    biological_replicates = unique(unlist(Seu[[biological_replicate_col]]))
    
    mat = matrix(0, nrow = nrow(Seu@assays$RNA@counts), ncol = length(cluster_u) * length(biological_replicates))
    rownames(mat) = rownames(Seu@assays$RNA@counts)
    
    n = 0
    names_mat =c()
    for(i in cluster_u){
        for(b in biological_replicates){
            n = n+1
            cells = colnames(Seu)[which(Seu$seurat_clusters == i & unlist(Seu[[biological_replicate_col]]) == b)]
            if(length(cells) > 25) {
                mat[,n] = Matrix::rowSums(Seu@assays$RNA@counts[,cells])
                names_mat = c(names_mat,paste0(i,"_",b))
            }
        }
    }
    mat = mat[,which(Matrix::colSums(mat) > 0)]
    colnames(mat) = names_mat
    return(mat)
}

#' Find Differentiated Clusters scRNA
#'
#' @details 
#'
#'
#' @param Seu A Seurat object containing scRNA dataset with 'seurat_clusters' 
#' column.
#' @param diff_method A character specifying the method to use for differential
#' analysis. Must be either "edgeR-LRT" or "Seurat".
#' @param biological_replicate_col Optional. If diff_method == "edgeR-LRT",
#' a character specifying a column of the Seurat object corresponding to the 
#' replicates or batches of the dataset in order to take in  account 
#' biological/technical noise for differential analysis. If NULL, will create 
#' 3 random layers of pseudo replicates.
#' @param qval.th A numeric specifying the adjusted p-value threshold above 
#' which a feature is defined as significantly differential.
#' @param logFC.th A numeric specifying the log Fold-Change of activation
#' (% cells activated) threshold above which a feature is defined as
#' significantly differential.
#' @param limit An integer specifying the minimum number of features required 
#' for a subcluster to be called 'true' subcluster.
#' @param cluster_of_origin A character specifying the name of the cluster of 
#' origin that will be concatenated before the name of true subclusters. 
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#'
#' @return
#' @export
#' 
#' 
#' @examples
find_differentiated_clusters_scRNA <- function(Seu,
                                               diff_method = c("edgeR-LRT", "Seurat")[1],
                                               logFC.th = log2(1.5),
                                               qval.th = 0.01,
                                               biological_replicate_col = NULL,
                                               min.pct = 0.1,
                                               limit = 5,
                                               cluster_of_origin = "Omega",
                                               min_frac_cell_assigned = 0.1
                                               ){
    
    # If Seurat is a SingleCellExperiment or else, try convert it to Seurat
    if(!is(Seu, "Seurat")) {
        Seu = Seurat::as.Seurat(Seu)
        if(is.null(Seu$seurat_clusters) & !is.null(Seu$seurat_clusters)) 
            Seu$seurat_clusters = Seu$cell_cluster else stop("Need 'seurat_clusters' or 'cell_cluster', column.")
        Seurat::Idents(Seu) = Seu$seurat_clusters
    }
    
    set.seed(47)
    
    # Starting list of clusters to re-cluster
    cluster_u = unique(Seu$seurat_clusters)
    diffmat_n = data.frame(n_differential = rep(0,length(cluster_u)),
                           cluster_of_origin = cluster_of_origin,
                           subcluster = cluster_u,
                           subcluster_after_diff = cluster_u)
    
    if(diff_method == "Seurat") {
        res = Seurat::FindAllMarkers(Seu, logfc.threshold = logFC.th, return.thresh = qval.th,
                                     test.use = "wilcox", min.pct = min.pct, only.pos = TRUE)
    } else {
        res = pseudobulk_LRT_edgeR(Seu, biological_replicate_col)
    }
    
    n_cell_assigned = 0
    for(i in seq_along(cluster_u)){
        group_cells = colnames(Seu)[which(Seu$seurat_clusters %in% cluster_u[i])]
        res. = res[which(res$cluster == cluster_u[i]),]
        res. = res.[which(res.$avg_log2FC > logFC.th &
                              res.$p_val_adj < qval.th &
                              res.$pct.1 > min.pct),]
        diffmat_n$n_differential[i] = nrow(res.)
        cat(cluster_u[i],"- Found",diffmat_n$n_differential[i], "differential regions.\n")
        if(diffmat_n$n_differential[i] < limit){
            cat(cluster_u[i], " cluster has less than", limit, "enriched features.\nAssigning the cells to cluster of origin.\n")
            diffmat_n$subcluster_after_diff[i] = cluster_of_origin
        } else{
            n_cell_assigned = n_cell_assigned + length(group_cells)
            diffmat_n$subcluster_after_diff[i] = paste0(cluster_of_origin, ":", cluster_u[i])
        }
    }
    passing = TRUE
    
    cat("Finished finding differences - ", n_cell_assigned/ ncol(Seu), " fraction of cells were assigned.\n")
    if((n_cell_assigned/ ncol(Seu)) < min_frac_cell_assigned){
        cat("Not enough cells were assigned - not clustering.\n")
        passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}


#'  Iterative Differential Clustering scRNA
#'
#' @param Seu A Seurat object preprocessed with Seurat.
#' @param output_dir The output directory in which to plot and save objects.
#' @param plotting A logical specifying wether to save the plots or not.
#' @param saving A logical specifying wether to save the data or not.
#' @param biological_replicate_col Optional. If diff_method == "edgeR-LRT",
#' a character specifying a column of the Seurat object corresponding to the 
#' replicates or batches of the dataset in order to take in  account 
#' biological/technical noise for differential analysis. If NULL, will create 
#' 3 random layers of pseudo replicates.
#' @param diff_method A character specifying the method to use for differential
#' analysis. Must be either "edgeR-LRT" or "Seurat".
#' @param logFC.th A numeric specifying the log Fold-Change of activation
#' (% cells activated) threshold above which a feature is defined as
#' significantly differential.
#' @param min.pct Minimum percentage of cells to be active in the cells of the
#' cluster to consider a feature as potentially significantly differential.
#' @param limit An integer specifying the minimum number of significantly 
#' enriched / depleted features required in order for a subcluster to be called
#' a 'true' subcluster
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param k An integer specifying the number of nearest neighbors to use for 
#' the Louvain clustering at each iteration.
#' @param resolution A numeric specifying the resolution to use for the Louvain
#' clustering at each iteration.
#' @param nPCA  An integer specifying the number of first PC to keep in the 
#' dimensionality reduction step
#'
#' @return A character vector containing the assignation of cells to clusters.
#' If saving is true, also saves list of differential analyses, differential 
#' analyses summaries and embeddings for each re-clustered cluster.
#'
#' @export
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom grDevices colors
#' @importFrom qs qsave
#' 
#' @examples
iterative_differential_clustering_scRNA <- function(
    Seu,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    biological_replicate_col = "Target",
    nPCA = 50,
    diff_method = c("edgeR-LRT", "Seurat")[1],
    logFC.th = log2(2),
    qval.th = 0.01,
    min.pct = 0.1,
    min_frac_cell_assigned = 0.25,
    limit = 10,
    k = 100,
    resolution = 0.5
){
    
    set.seed(47)
    
    # Colors for the plot
    if (is.null(color)){
        color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
        color <- sample(color, 400)
    } else {
        stopifnot(is.character(color))
        if(length(color) < 20) warning("IDclust::iterative_differential_clustering_scRNA - ",
                                       "The color vector might be too short.",
                                       "Please precise more colors.")
    }
    
    # For the first partition, try to find very low level clusters (e.g. low
    # resolution, high number of neighbors)
    Seu = Seurat::FindNeighbors(Seu, reduction = "pca",  k.param = 50, dims = 1:nPCA, verbose = F)
    Seu = Seurat::FindClusters(Seu, algorithm = 2, resolution = starting.resolution, random.seed = 47, verbose = FALSE)
    
    Seu$seurat_clusters <- paste0("A",as.numeric(Seu$seurat_clusters))
    Seurat::Idents(Seu) = Seu$seurat_clusters
    
    # Starting list of clusters to re-cluster
    cluster_u = unique(Seu$seurat_clusters)
    
    # Find initial differences (even if there are none, the initial clusters
    # are always considered as true clusters).
    cat("Initial round of clustering - limit of differential genes set to 0",
        " for this first round only.\n")
    DA = find_differentiated_clusters_scRNA(
        Seu, logFC.th = logFC.th, qval.th = qval.th,
        min.pct = min.pct, diff_method = diff_method,
        biological_replicate_col = biological_replicate_col,
        min_frac_cell_assigned = min_frac_cell_assigned, 
        cluster_of_origin =  "Omega",
        limit = 0)
    
    diffmat_n = DA$diffmat_n
    
    # List of embeddings
    list_embeddings = list(Seu@reductions$pca@cell.embeddings )
    names(list_embeddings)[1] = paste(unique(Seu$seurat_clusters), collapse = "_")
    
    # List of differential analyses
    list_res = list(DA$res)
    names(list_res)[1] = "Omega"
    
    # List of DA summaries
    list_diffmat =  list(diffmat_n)
    names(list_diffmat)[1] = "Omega"
    
    iteration = 0
    gc()
    
    while(iteration < nrow(differential_summary_df)){
        
        if(plotting == TRUE){
            # Plot initial
            png(file.path(output_dir, paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
            print(
                Seurat::DimPlot(Seu, group.by =  "seurat_clusters", cols = color) 
            )
            dev.off()
        }
        iteration = iteration + 1
        cat("Doing iteration number ", iteration,"...\n")
        
        
        partition_cluster_of_origin = differential_summary_df$cluster[iteration]
        if(!(partition_cluster_of_origin  %in% differential_summary_df$cluster[duplicated(differential_summary_df$cluster)])){
            
            partition_depth = which(LETTERS == substr(gsub(".*:","",partition_cluster_of_origin),1,1)) + 1
            
            # Select first cluster
            Seu. = Seu[, which(Seu$seurat_clusters %in%  partition_cluster_of_origin)]
            cat("Re-calculating PCA and subclustering for cluster", partition_cluster_of_origin,".\n")
            
            if(ncol(Seu.) > 100){
                # Rerun PCA
                Seu. = Seurat::FindVariableFeatures(Seu., selection.method = "vst", verbose = FALSE)
                Seu. = Seurat::ScaleData(Seu., verbose = FALSE)
                Seu. = Seurat::RunPCA(Seu., features = Seurat::VariableFeatures(object = Seu.), verbose = FALSE)
                Seu. = Seurat::FindNeighbors(Seu., use.dimred = "pca", verbose = FALSE)
                Seu. = Seurat::FindClusters(Seu., algorithm = 2, resolution = 0.1,  random.seed = 47, verbose = FALSE)
                Seu.$seurat_clusters = paste0(LETTERS[partition_depth], as.numeric(Seu.$seurat_clusters))
                Seurat::Idents(Seu.) = Seu.$seurat_clusters
                
                if(plotting == TRUE){
                    png(file.path(output_dir, paste0(partition_cluster_of_origin,"_raw.png")), width = 1400, height = 1200, res = 200)
                    print(
                        Seurat::DimPlot(Seu., reduction = "pca") + ggtitle(paste0(partition_cluster_of_origin, " - raw"))
                    )
                    dev.off()
                }
                clusters = Seu.$seurat_clusters 
                cluster_u = unique(clusters)
                
                
                if(length(cluster_u) > 1 ){
                    cat("Found", length(cluster_u),"subclusters.\n")
                    if(verbose) print(table(clusters))
                    
                    ## Differential analysis
                    DA = find_differentiated_clusters_scRNA(Seu., diff_method = diff_method,
                                                            logFC.th = logFC.th, qval.th = qval.th,
                                                            min.pct = min.pct,  limit = limit,
                                                            biological_replicate_col = biological_replicate_col,
                                                            min_frac_cell_assigned = min_frac_cell_assigned, 
                                                            cluster_of_origin =  partition_cluster_of_origin)
                    
                    # Retrieve DA results
                    diffmat_n = DA$diffmat_n
                    list_res[[partition_cluster_of_origin]] = DA$res
                    list_embeddings[[partition_cluster_of_origin]] = Seu.@reductions$pca@cell.embeddings
                    list_diffmat[[partition_cluster_of_origin]] = diffmat_n
                    
                    # If more than 'min_frac_cell_assigned' of the cells were assigned
                    # to 'true' subclusters (with marker features)
                    if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                        
                        # Add the new sublclusters to the list of clusters
                        df = data.frame("cluster" = diffmat_n$subcluster_after_diff, n_differential = diffmat_n$n_differential)
                        differential_summary_df = rbind(differential_summary_df, df)
                        Seu.$seurat_clusters = diffmat_n$subcluster_after_diff[match(Seu.$seurat_clusters, diffmat_n$subcluster)]
                        Seu$seurat_clusters[match(colnames(Seu.), colnames(Seu))] = Seu.$seurat_clusters
                        
                        Seurat::Idents(Seu.) = Seu.$seurat_clusters
                        
                        if(plotting == TRUE){
                            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                            print(
                                Seurat::DimPlot(Seu., reduction = "pca") + ggtitle(paste0(partition_cluster_of_origin, " - true"))
                            )
                            dev.off()
                        }
                    } else{
                        cat("Found 2 subcluster for", partition_cluster_of_origin," but with no difference. Maximum clustering reached...\n")
                    }
                    
                } else{
                    cat("Found only 1 subcluster for", partition_cluster_of_origin,". Maximum clustering reached...\n")
                }
            }
            gc()
            
        } else{
            cat(partition_cluster_of_origin,"- This cluster was formed by 'unassigned' cells, not clustering it further...\n")
        }
    }
    
    cat("\n\n\n##########################################################\nFinished !\nFound a total of", length(unique(Seu$seurat_clusters)),"clusters after",iteration ,"iterations.",
        "\nThe average cluster size is ",floor(mean(table(Seu$seurat_clusters)))," and the median is",floor(median(table(Seu$seurat_clusters))),".",
        "\nThe number of initital clusters not subclustered is ",length(grep(":", unique(Seu$seurat_clusters),invert = TRUE)),".",
        "\n##########################################################\n")
    
    
    ## Saving results
    if(saving){
        # List of differential features for each re-clustering
        qs::qsave(list_res, file.path(output_dir, "IDC_DA.qs"))
        
        # List of embedding of each re-clustered cluster
        qs::qsave(list_embeddings, file.path(output_dir, "IDC_embeddings.qs"))
        
        # List of summaries of the number of differential features for each re-clustering
        qs::qsave(list_diffmat, file.path(output_dir, "IDC_summaries_DA.qs"))
        
        # Final SingleCellExperiment with the clusters found by IDC 
        qs::qsave(Seu, file.path(output_dir, "Seu_IDC.qs"))
    }
    
    return(Seu$seurat_clusters)
}

