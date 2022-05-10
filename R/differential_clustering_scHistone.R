
#' Find Differentiated Clusters scEpigenomics
#'
#' @details  Find significantly differential features between the given set
#' of clusters (within the 'cell_cluster' column of the SingleCellExperiment).
#' For each cluster, if enough differences are found, mark the cluster as a 
#' 'true' subcluster and gives it the alias 'cluster_of_origin:cluster'. 
#' The function will use by default ChromSCape::differential_activation function
#' to define differential features.
#' 
#' @param scExp A SingleCellExperiment with 'cell_cluster' column filled with 
#' cluster assignations.  
#' @param qval.th A numeric specifying the adjusted p-value threshold above 
#' which a feature is defined as significantly differential.
#' @param min.frac A numeric between 0 and 1 specifying the minimum fraction of
#' cells activated for a given feature to consider it as differential.
#' @param logFC.th A numeric specifying the log Fold-Change of activation
#' (% cells activated) threshold above which a feature is defined as
#' significantly differential.
#' @param limit An integer specifying the minimum number of features required 
#' for a subcluster to be called 'true' subcluster.
#' @param limit_by_proportion Optional. A data.frame containing 3 columns - 
#' ncells, mean_n_differential, sd_n_differential - that reflect the number of 
#' false positive expected for a given cluster size.
#'  (See  \link{calculate_FDR_scEpigenomics}).
#' @param cluster_of_origin A character specifying the name of the cluster of 
#' origin that will be concatenated before the name of true subclusters. 
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#'
#' @return A list containing :
#'  * "diffmat_n" - A data.frame containing the number of differential regions
#'  foun per cluster and the new assignations of the subclusters.
#'  * "res" - A data.frame containing the differential analysis.
#'  * "passing_min_pct_cell_assigned" - A boolean indicating if enough cells
#'   were  assigned
#' @export
#' 
#' @import dplyr
#' 
#' @examples
#'
find_differentiated_clusters_scEpigenomics <- function(scExp,
                                                       qval.th = 0.1,
                                                       min.frac = 0.1,
                                                       logFC.th = 0.5,
                                                       limit = 5,
                                                       limit_by_proportion = NULL,
                                                       min_frac_cell_assigned = 0.1,
                                                       cluster_of_origin = "Omega"){
    cluster_u = names(sort(table(scExp$cell_cluster)))
    diffmat_n = data.frame(n_differential = rep(0,length(cluster_u)),
                           cluster_of_origin = cluster_of_origin,
                           subcluster = cluster_u,
                           subcluster_after_diff = cluster_u)
    
    res = ChromSCape::differential_activation(scExp = scExp, group_by = "cell_cluster", verbose = TRUE)
    gc()
    n_cell_assigned = 0
    for(i in 1:(length(cluster_u))){
        group = cluster_u[i]
        group_cells = colnames(scExp)[scExp$cell_cluster == group]
        diffmat_n$n_differential[i] = length(which(
            res[,paste0("logFC.", group)] > logFC.th  & 
                res[,paste0("qval.", group)] < qval.th  & 
                res[,paste0("group_activation.", group)] > min.frac
        ))
        
        cat(cluster_u[i]," - found", diffmat_n$n_differential[i], "enriched features.\n")
        if(is.null(limit_by_proportion)){
            if(diffmat_n$n_differential[i] < limit){
                cat(cluster_u[i], " cluster has less than", limit, "enriched features.\nAssigning the cells to cluster of origin.\n")
                diffmat_n$subcluster_after_diff[i] = cluster_of_origin
            } else{
                n_cell_assigned = n_cell_assigned + length(group_cells)
                diffmat_n$subcluster_after_diff[i] = paste0(cluster_of_origin, ":", cluster_u[i])
            }
        } else{
            index = which.min(abs(limit_by_proportion$ncells - length(group_cells)))
            limit. = max(limit, 
                         limit_by_proportion$mean_n_differential[index][1] + 2 * limit_by_proportion$sd_n_differential[index][1])
            if(diffmat_n$n_differential[i] < limit.){
                cat(cluster_u[i], " cluster has less than", limit., "enriched features.\nAssigning the cells to cluster of origin.\n")
                diffmat_n$subcluster_after_diff[i] = cluster_of_origin
            } else{
                n_cell_assigned = n_cell_assigned + length(group_cells)
                diffmat_n$subcluster_after_diff[i] = paste0(cluster_of_origin, ":", cluster_u[i])
            }
        }
    }
    passing = TRUE
    
    cat("Finished finding differences - ", n_cell_assigned/ ncol(scExp), " fraction of cells were assigned.\n")
    if((n_cell_assigned/ ncol(scExp)) < min_frac_cell_assigned){
        cat("Not enough cells were assigned - not clustering.\n")
        passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}


#' Title
#'
#' @param scExp A SingleCellExperiment object 
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
#' @return
#' @export
#'
#' @examples
calculate_FDR_scEpigenomics <- function(scExp,
                                        limit = 5,
                                        qval.th = 0.01,
                                        logFC.th = log2(2),
                                        min.pct = 0.01,
                                        min_frac_cell_assigned = 0.1,
                                        cluster_of_origin = "Omega",
                                        nThreads = 10){
    if(!requireNamespace("doParallel")){
        stop("IDclust::calculate_FDR_scEpigenomics - Please install ",
             "doParallel to run FDR in parallel")
    }
    if(!requireNamespace("foreach")){
        stop("IDclust::calculate_FDR_scEpigenomics - Please install ",
             "foreach to run FDR in parallel")
    }
    
    myCluster <- parallel::makeCluster(nThreads, type = "FORK")
    doParallel::registerDoParallel(myCluster)
    FDR_iterations = foreach::foreach(i = 1:10) %dopar% {
        set.seed(i * 47)
        scExp$cell_cluster = sample(scExp$cell_cluster, length(scExp$cell_cluster), replace = F)
        DA <- run_real_one_vs_all_comparisons_activation_scExp(scExp = scExp,
                                                               qval.th = qval.th,
                                                               min.pct = min.pct.activation,
                                                               logFC.th = log2(FC.th),
                                                               min_frac_cell_assigned = min_frac_cell_assigned,
                                                               limit = limit,
                                                               cluster_of_origin = partition_cluster_of_origin)
        df = DA$diffmat_n[,-4]
        df$iteration = i
        df
    }  

    FDR = do.call("rbind", FDR_iterations)
    
    parallel::stopCluster(myCluster)
    gc()
    
    return(FDR)
}

#' Iterative Differential Clustering scEpigenomics
#'
#' @description Main function of the IDclust package. Provided a 
#' SingleCellExperiment pre-processed with ChromSCape, will find biologically 
#' relevant clusters by iteratively re-clustering and re-processing clusters. 
#' At each iteration, subclusters having enough significantly enriched features 
#' compared to other subclusters are defined as 'true' subclusters. Others are 
#' assigned to parent clusters. The algorithm will stop when no more 'true' 
#' subclusters are found. 
#' 
#' This method ensure that each cluster found in this unsupervised way have
#' significant biological differences, based on the user defined thresholds.
#' 
#' @details The default differential analysis used is the 
#' differential_activation function from ChromSCape package. This function 
#' compares the % of active cells in the cluster versus the rest of cells and 
#' perform a Chi-squared test to calculate p-values.  
#' 
#'
#' @param scExp A SingleCellExperiment object 
#' @param output_dir The output directory in which to plot and save objects.
#' @param plotting A logical specifying wether to save the plots or not.
#' @param saving A logical specifying wether to save the data or not.
#' @param limit An integer specifying the minimum number of significantly 
#' enriched / depleted features required in order for a subcluster to be called
#' a 'true' subcluster
#' @param FC.th  A numeric specifying the fold change of activation above/below
#' which a feature is considered as significantly differential.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential.
#' @param quantile.activation A numeric between 0 and 1 specifying the quantile
#' of global activation to take as minimal percentage of activation for the 
#' differential analysis. Increasing this value will decrease the number of 
#' differential features.
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param k An integer specifying the number of nearest neighbors to use for 
#' the Louvain clustering at each iteration.
#' @param resolution A numeric specifying the resolution to use for the Louvain
#' clustering at each iteration.
#' @param starting.resolution A numeric specifying the resolution to use for the 
#' Louvain clustering of the first iteration. It is recommended to set it quite
#' low in order to have few starting clusters.
#' @param limit_by_proportion Optional. A data.frame containing 3 columns - 
#' ncells, mean_n_differential, sd_n_differential - that reflect the number of 
#' false positive expected for a given cluster size.
#' @param frac_feature A numeric between 0 and 1 specifying the fraction of 
#' features to keep at each iteration when re-processing the subset of cells.
#' @param nPCA  An integer specifying the number of first PC to keep in the 
#' dimensionality reduction step.
#' @param runFDR A logical indicating wether to also run the FDR calculation at
#' each iteration. This operation is optional and requires time x CPU x RAM.
#' @param nThreads If runFDR==TRUE. An integer specifying of threads to use
#' for the calculation of the FDR.
#' @param color Set of colors to use for the coloring of the clusters. This must
#' contains enough colors for each cluster (minimum 20 colors, but 100 colors
#' at least is recommended, based on the dataset).
#' @param verbose A logical specifying wether to print.
#'
#'
#' @return A character vector containing the assignation of cells to clusters.
#' If saving is true, also saves list of differential analyses, differential 
#' analyses summaries and embeddings for each re-clustered cluster. If runFDR is 
#' TRUE, also saves the list of FDR for each re-clusterd cluster.
#' 
#' @export
#' 
#' @importFrom  Matrix Matrix rowSums
#' @importFrom qs qsave
#' @importFrom  SummarizedExperiment assays 
#' @importFrom  grDevices colors 
#' @importFrom  SingleCellExperiment reducedDim  
#' @import ggplot2
#' @import dplyr
#' 
#' @examples
#' 
iterative_differential_clustering_scEpigenomics <- function(
    scExp,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    nPCA = 10,
    frac_feature = 1,
    quantile.activation = 0.7,
    min_frac_cell_assigned = 0.1,
    FC.th = 2,
    qval.th = 0.1,
    limit = 5,
    k = 100,
    starting.resolution = 0.1,
    resolution = 0.8,
    runFDR = FALSE,
    limit_by_proportion = NULL,
    color = NULL,
    nThreads = 10,
    verbose = TRUE
){
    set.seed(47)
    
    # For the first partition, try to find very low level clusters (e.g. low
    # resolution, high number of neighbors)
    scExp = ChromSCape::find_clusters_louvain_scExp(scExp,
                                                    k = 100,
                                                    resolution = starting.resolution,
                                                    use.dimred = "PCA")
    scExp$cell_cluster = gsub("C","A", scExp$cell_cluster)
    
    # Calculate the average % cells activated in a feature
    # and return a level of activation based on a given decile
    binmat = Matrix::Matrix((SummarizedExperiment::assays(scExp)$counts > 0) + 0, sparse = TRUE)
    min.frac.activation = quantile(Matrix::rowSums(binmat) / ncol(binmat), quantile.activation)
    
    
    # Find initial differences (even if there are none, the initial clusters
    # are always considered as true clusters).
    cat("Initial round of clustering - limit of differential genes set to 0",
        " for this first round only.\n")
    DA = find_differentiated_clusters_scEpigenomics(scExp = scExp,
                                                    qval.th = qval.th,
                                                    min.frac = min.frac.activation,
                                                    logFC.th = log2(FC.th),
                                                    min_frac_cell_assigned = min_frac_cell_assigned,
                                                    limit = 0,
                                                    limit_by_proportion = NULL,
                                                    cluster_of_origin = "Omega")
    
    # Starting list of clusters to re-cluster
    differential_summary_df = data.frame("cluster" = unique(scExp$cell_cluster),
                                         "n_differential" = 10) # we assume the initial clusters are differential
    
    # Colors for the plot
    if (is.null(color)){
        color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
        color <- sample(color, 400)
    } else {
        stopifnot(is.character(color))
        if(length(color) < 20) warning("IDclust::iterative_differential_clustering - ",
                                       "The color vector might be too short.",
                                       "Please precise more colors.")
    }
    scExp$cell_cluster_color = NULL
    
    # List of embeddings
    list_embeddings = list(SingleCellExperiment::reducedDim(scExp,"PCA"))
    names(list_embeddings)[1] = paste(unique(scExp$cell_cluster), collapse = "_")
    
    # List of marker features
    list_res = list(DA$res)
    names(list_res)[1] = "Omega"
    
    # List of summary of marker features
    list_diffmat = list(DA$diffmat_n)
    names(list_diffmat)[1] = "Omega"
    
    # List of FDR of each clusters
    FDR_list = list()
    
    iteration = 0
    gc()
    
    
    # Run IDC until no more clusters can be re-clustered into meaningful clusters
    while(iteration < nrow(differential_summary_df)){
        
        if(plotting == TRUE){
            # Plot each iteration of the algorithm
            png(file.path(output_dir, paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
            print(
                ChromSCape::plot_reduced_dim_scExp(scExp, reduced_dim = "UMAP", color_by = "cell_cluster",
                                                   downsample = 50000, size = 0.35, transparency = 0.75) +
                    ggtitle(paste0("Initital Clustering")) + theme(legend.position = "none")
            )
            dev.off()
        }
        
        iteration = iteration + 1
        
        cat("Doing iteration number ", iteration,"...\n")
        
        # Name of the cluster of origin
        partition_cluster_of_origin = differential_summary_df$cluster[iteration]
        
        # Make sure that we do not re-cluster a 'parent cluster',
        # e.g. cells that were unassigned and therefore assigned to the parent
        # cluster
        if(!(partition_cluster_of_origin  %in% differential_summary_df$cluster[duplicated(differential_summary_df$cluster)])){
            
            # Letter assigned to the 'partition depth' (e.g. A, B, C...)
            partition_depth = which(LETTERS == substr(gsub(".*:","",partition_cluster_of_origin),1,1)) + 1
            
            # Select only cells from the given cluster
            scExp. = scExp[, which(scExp$cell_cluster %in%  partition_cluster_of_origin)]
            
            cat("Re-calculating PCA and subclustering for cluster", partition_cluster_of_origin,
                "with",ncol(scExp.),"cells.\n")
            
            # Re-cluster a cluster only if there are more than 100 cells
            if(ncol(scExp.) > 100){
                
                # Re-run TFIDF and PCA and clusters using Louvain algorithm 
                scExp. = ChromSCape::find_top_features(scExp.,n = floor(frac_feature * nrow(scExp.)),
                                                       keep_others = FALSE)
                scExp. = ChromSCape::normalize_scExp(scExp., type = "TFIDF")
                scExp. = ChromSCape::reduce_dims_scExp(scExp., dimension_reductions = "PCA",
                                                       n = nPCA, remove_PC = "Component_1")
                scExp. = ChromSCape::find_clusters_louvain_scExp(scExp., k = 50, resolution =  resolution,
                                                                 use.dimred = "PCA")
                
                # Give the cluster the letter corresponding to the partition depth
                scExp.$cell_cluster <- paste0(LETTERS[partition_depth],gsub("C", "", scExp.$cell_cluster))
                
                # Plot the 'raw' PCA on the given cluster before re-assigning cells
                if(plotting == TRUE){
                    png(file.path(output_dir, paste0(partition_cluster_of_origin,"_raw.png")), width = 1400, height = 1200, res = 200)
                    print(
                        ChromSCape::plot_reduced_dim_scExp(scExp., color_by = "cell_cluster", reduced_dim = "PCA") + 
                            ggtitle(paste0(partition_cluster_of_origin, " - raw"))
                    )
                    dev.off()
                }
                clusters = scExp.$cell_cluster 
                cluster_u = unique(clusters)
                
                if(length(cluster_u) > 1 ){
                    cat("Found", length(cluster_u),"subclusters.\n")
                    if(verbose) print(table(clusters))
                    
                    # Find differentiated clusters
                    DA = find_differentiated_clusters_scEpigenomics(scExp = scExp.,
                                                                    qval.th = qval.th,
                                                                    min.frac = min.frac.activation,
                                                                    logFC.th = log2(FC.th),
                                                                    min_frac_cell_assigned = min_frac_cell_assigned,
                                                                    limit = limit,
                                                                    limit_by_proportion = limit_by_proportion,
                                                                    cluster_of_origin = partition_cluster_of_origin)
                    gc()
                    
                    # Retrieve DA results
                    diffmat_n = DA$diffmat_n
                    list_res[[partition_cluster_of_origin]] = DA$res
                    list_embeddings[[partition_cluster_of_origin]] = SingleCellExperiment::reducedDim(scExp., "PCA")
                    list_diffmat[[partition_cluster_of_origin]] = diffmat_n
                    
                    # If more than 'min_frac_cell_assigned' of the cells were assigned
                    # to 'true' subclusters (with marker features)
                    if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                        
                        # Add the new sublclusters to the list of clusters
                        df = data.frame("cluster" = diffmat_n$subcluster_after_diff, n_differential = diffmat_n$n_differential)
                        differential_summary_df = rbind(differential_summary_df, df)
                        scExp.$cell_cluster = diffmat_n$subcluster_after_diff[match(scExp.$cell_cluster, diffmat_n$subcluster)]
                        scExp$cell_cluster[match(colnames(scExp.), colnames(scExp))] = scExp.$cell_cluster
                        
                        # Calculate FDR by mixing the cluster assignation 
                        #  10 times and calculate the number of false positive
                        if(runFDR == TRUE){
                            
                            FDR = calculate_FDR_scEpigenomics(scExp,
                                                              limit = limit,
                                                              qval.th = qval.th,
                                                              logFC.th = log2(FC.th),
                                                              min.pct = min.pct.activation,
                                                              min_frac_cell_assigned = min_frac_cell_assigned,
                                                              cluster_of_origin = "Omega",
                                                              nThreads = nThreads)
                            
                            FDR_list[[partition_cluster_of_origin]] = FDR
                            FDR_val = length(which(FDR_list[[partition_cluster_of_origin]]$n_differential >= 5)) / nrow(FDR_list[[partition_cluster_of_origin]])
                            gc()
                        }
                        
                        tit = ifelse(runFDR, paste0(" - FDR - ", FDR_val), "")
                        if(plotting == TRUE){
                            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                            print(
                                plot_reduced_dim_scExp(scExp., color_by = "cell_cluster", reduced_dim = "PCA") + 
                                    ggtitle(paste0(partition_cluster_of_origin, tit))
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
    
    cat("\n\n\n##########################################################\nFinished !\nFound a total of", length(unique(scExp$cell_cluster)),"clusters after",iteration ,"iterations.",
        "\nThe average cluster size is ", floor(mean(table(scExp$cell_cluster)))," and the median is",floor(median(table(scExp$cell_cluster))),".",
        "\nThe number of initital clusters not subclustered is ",length(grep(":", unique(scExp$cell_cluster),invert = TRUE)),".",
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
        qs::qsave(scExp, file.path(output_dir, "scExp_IDC.qs"))
    }
    # List of FDR for each re-clustering
    if(runFDR) qs::qsave(FDR_list, file.path(output_dir, "FDR_list.qs"))
    
    return(scExp$cell_cluster)
}

