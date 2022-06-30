#' @importFrom grDevices dev.off png
#' @importFrom methods is 
#' @importFrom stats median model.matrix quantile
NULL


#' Find Differentiated Clusters
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @rdname find_differentiated_clusters
#' @export find_differentiated_clusters
#'
#' @examples
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' # Find differentiated clusters in Seurat object using 
#' # edgeR_pseudobulk_LRT function
#' data("Seu")
#'
#' DA = find_differentiated_clusters(
#'     Seu,
#'     differential_function = differential_edgeR_pseudobulk_LRT,
#'     logFC.th = log2(1.5),
#'     qval.th = 0.01,
#'     by = "seurat_clusters",
#'     limit = 5,
#'     cluster_of_origin = "Omega",
#'     min_frac_cell_assigned = 0.1,
#'     verbose = TRUE
#' )
#' 
#' # Summary of differential genes per cluster
#' head(DA$diffmat_n)
#' 
#' # Differential analysis
#' head(DA$res)
#' 
#' # Did the clustering pass the minimum percent of cell assigned threshold ?
#' print(DA$passing_min_pct_cell_assigned)
#' 
#' # Find differentiated clusters in Seurat object using Seurat function
#' data("Seu")
#' 
#' DA = find_differentiated_clusters(
#'     Seu,
#'     differential_function = differential_Seurat,
#'     logFC.th = log2(1.5),
#'     qval.th = 0.01,
#'     by = "seurat_clusters",
#'     limit = 5,
#'     cluster_of_origin = "Omega",
#'     min_frac_cell_assigned = 0.1,
#'     verbose = TRUE
#' )
#' 
#' # Find differentiated clusters in Seurat object using Seurat function, 
#' # passing additional arguments to differential_Seurat and thus 
#' # Seurat::FindAllMarkers funtion.
#' data("Seu")
#' 
#' DA = find_differentiated_clusters(
#'     Seu,
#'     differential_function = differential_Seurat,
#'     logFC.th = log2(1.5),
#'     qval.th = 0.01,
#'     by = "seurat_clusters",
#'     limit = 5,
#'     cluster_of_origin = "Omega",
#'     min_frac_cell_assigned = 0.1,
#'     verbose = TRUE,
#'     test.use = "roc" # additional argument
#' )
#' }
#' # Find differentiated clusters in SingleCellExperiment object using 
#' # differential_ChromSCape function.
#' if(requireNamespace("ChromSCape", quietly=TRUE)){
#' data("scExp")
#' scExp = ChromSCape::find_clusters_louvain_scExp(scExp, 
#' resolution = 0.1)
#' DA = find_differentiated_clusters(
#'     scExp,
#'     differential_function = differential_ChromSCape,
#'     logFC.th = log2(5),
#'     qval.th = 0.01,
#'     by = "IDcluster",
#'     limit = 5,
#'     cluster_of_origin = "Omega",
#'     min_frac_cell_assigned = 0.1,
#'     verbose = TRUE,
#' )
#' }
find_differentiated_clusters <- function(object, ...) {
    UseMethod(generic = 'find_differentiated_clusters', object = object)
    
}

#' iterative_differential_clustering
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @rdname iterative_differential_clustering
#' @export iterative_differential_clustering
#'
#' @examples
#' # Clustering of Seurat scRNA object (Paired-Tag)
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' 
#' data("Seu", package = "IDclust")
#' set.seed(47)
#' Seu = iterative_differential_clustering(Seu, saving = FALSE, plotting =FALSE,
#' logFC.th = 0.2, qval.th = 0.1)
#' 
#' }
#' 
#' # Clustering of scExp scH3K27ac object (Paired-Tag)
#' if(requireNamespace("ChromSCape", quietly=TRUE)){
#' 
#' data("scExp", package = "IDclust")
#' set.seed(47)
#' scExp = iterative_differential_clustering(scExp, saving = FALSE, plotting =FALSE,
#' logFC.th = 0.5, qval.th = 0.01)
#' 
#' }
iterative_differential_clustering <- function(object, ...) {
    UseMethod(generic = 'iterative_differential_clustering', object = object)
    
}

#' @details  Find significantly differential features between the given set
#' of clusters (within the 'IDcluster' column of the SingleCellExperiment).
#' For each cluster, if enough differences are found, mark the cluster as a 
#' 'true' subcluster and gives it the alias 'cluster_of_origin:cluster'. 
#' The function will use by default [ChromSCape::differential_activation()]
#' function to define differential features.
#' 
#' @param object A SingleCellExperiment with 'IDcluster' column filled with 
#' cluster assignations.  
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster.
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential 
#' passed to the differential_function.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential passed to the
#' differential_function.
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param limit An integer specifying the minimum number of features required 
#' for a subcluster to be called 'true' subcluster.
#' @param FP_linear_model Optional. A linear model (see [stats::lm()]) of 
#' the number of false positive expected for a given cluster size. The [lm_list]
#' list of linear models present in this package gives default values accross
#' multiple binsizes. (See  [calculate_FDR_scEpigenomics]).
#' @param cluster_of_origin A character specifying the name of the cluster of 
#' origin that will be concatenated before the name of true subclusters. 
#' @param min_cluster_size An integer specifying the minimum number of cells
#' in a cluster to consider it as a 'true' subcluster.
#' @param verbose A logical specifying wether to print.
#' @param ... Additional parameters passed to the differential_function. See 
#' [differential_ChromSCape()] for more information on additional
#' parameters for the default function.
#'
#' @seealso [differential_ChromSCape()]
#' 
#' @return A list containing :
#'  * "diffmat_n" - A data.frame containing the number of differential regions
#'  foun per cluster and the new assignations of the subclusters.
#'  * "res" - A data.frame containing the differential analysis.
#'  * "passing_min_pct_cell_assigned" - A boolean indicating if enough cells
#'   were  assigned
#' 
#' @import dplyr
#' @export 
#' @rdname find_differentiated_clusters
#' @exportS3Method  find_differentiated_clusters default
find_differentiated_clusters.default <- function(
    object,
    differential_function = differential_ChromSCape,
    by = "IDcluster",
    logFC.th = log2(1.5),
    qval.th = 0.01,
    min_frac_cell_assigned = 0.1,
    limit = 5,
    FP_linear_model = NULL,
    cluster_of_origin = "Omega",
    min_cluster_size = 30,
    verbose = TRUE,
    ...){
    cluster_u = names(sort(table(object[[by]])))
    diffmat_n = data.frame(n_differential = rep(0,length(cluster_u)),
                           cluster_of_origin = cluster_of_origin,
                           subcluster = cluster_u,
                           true_subcluster = cluster_u)
    
    res = differential_function(object, by = by, logFC.th = logFC.th, qval.th = qval.th, ...)
    gc()
    n_cell_assigned = 0
    for(i in 1:(length(cluster_u))){
        cluster = cluster_u[i]
        group_cells = colnames(object)[object[[by]] == cluster]
        if(length(group_cells)  >= min_cluster_size){
            
            res. = res[which(res$cluster == cluster),]
            
            diffmat_n$n_differential[i] = nrow(res.) 
            
            if(verbose) cat(cluster," - found", diffmat_n$n_differential[i], "enriched features.\n")
            
            if(is.null(FP_linear_model)){
                if(diffmat_n$n_differential[i] < limit){
                    if(verbose) cat(cluster, " cluster has less than", limit, "enriched features.\nAssigning the cells to cluster of origin.\n")
                    diffmat_n$true_subcluster[i] = cluster_of_origin
                } else{
                    n_cell_assigned = n_cell_assigned + length(group_cells)
                    diffmat_n$true_subcluster[i] = paste0(cluster_of_origin, ":", cluster)
                }
            } else{
              cluster_size <- data.frame(ncells=length(group_cells))
              conf = stats::predict(FP_linear_model, newdata = cluster_size, interval = 'confidence')
              confidence_limit = conf[1,"upr"]
              cat("Ncells = ",length(group_cells), " - Confidence = ", confidence_limit, "\n")
                limit. = max(limit,  confidence_limit)
                if(diffmat_n$n_differential[i] < limit.){
                    if(verbose) cat(cluster, " cluster has less than", limit., "enriched features.\nAssigning the cells to cluster of origin.\n")
                    diffmat_n$true_subcluster[i] = cluster_of_origin
                } else{
                    n_cell_assigned = n_cell_assigned + length(group_cells)
                    diffmat_n$true_subcluster[i] = paste0(cluster_of_origin, ":", cluster)
                }
            }
        } else{
            if(verbose) cat(cluster, " cluster has less than", min_cluster_size, " cells. \nAssigning the cells to cluster of origin.\n")
            diffmat_n$true_subcluster[i] = cluster_of_origin
        }
    }
    passing = TRUE
    
    if(verbose) cat("Finished finding differences - ", n_cell_assigned/ ncol(object), " fraction of cells were assigned.\n")
    if((n_cell_assigned/ ncol(object)) < min_frac_cell_assigned){
        if(verbose) cat("Not enough cells were assigned - not clustering.\n")
        passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}


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
#' [ChromSCape::differential_activation()] function. This function 
#' compares the % of active cells in the cluster versus the rest of cells and 
#' perform a Chi-squared test to calculate p-values.   
#' 
#'
#' @param object A SingleCellExperiment object.
#' @param output_dir The output directory in which to plot and save objects.
#' @param plotting A logical specifying wether to save the plots or not.
#' @param saving A logical specifying wether to save the data or not.
#' @param limit An integer specifying the minimum number of significantly 
#' enriched / depleted features required in order for a subcluster to be called
#' a 'true' subcluster
#' @param dim_red The name of the slot to save the dimensionality reduction at
#' each step.
#' @param vizualization_dim_red The name of the slot used for plotting. Must be
#' a valid slot present in \code{reducedDimNames(object)}.
#' @param processing_function A function that re-process the subset of clusters
#' at each step. It msut take in entry a  SingleCellExperiment object, \code{dim_red}
#'  and \code{n_dims} parameters and returns a SingleCellExperiment containing a cell
#' embedding. See [processing_ChromSCape] for the default function.
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential 
#' passed to the differential_function.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential passed to the
#' differential_function.
#' @param quantile.activation A numeric between 0 and 1 specifying the quantile
#' of global activation to take as minimal percentage of activation for the 
#' differential analysis. Increasing this value will decrease the number of 
#' differential features.
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param starting.resolution A numeric specifying the resolution to use for the 
#' Louvain clustering of the first iteration. It is recommended to set it quite
#' low in order to have few starting clusters.
#' @param starting.k An integer specifying the number of nearest neighbors to
#'  use for  the Louvain clustering of the first iteration. It is recommended 
#'  to set it quite high in order to have few starting clusters
#' @param resolution A numeric specifying the resolution to use for the Louvain
#' clustering at each iteration.
#' @param k An integer specifying the number of nearest neighbors to use for 
#' the Louvain clustering at each iteration.
#' @param FP_linear_model Optional. A linear model (see [stats::lm()]) of 
#' the number of false positive expected for a given cluster size. The [lm_list]
#' list of linear models present in this package gives default values accross
#' multiple binsizes. (See  [calculate_FDR_scEpigenomics]).
#' @param n_dims An integer specifying the number of first dimensions to keep 
#' in the dimensionality reduction step.
#' @param nThreads  An integer specifying of threads to use
#' for the calculation of the FDR.
#' @param color Set of colors to use for the coloring of the clusters. This must
#' contains enough colors for each cluster (minimum 20 colors, but 100 colors
#' at least is recommended, based on the dataset).
#' @param swapExperiment A character specifying an alternative experiment
#' (see [SingleCellExperiment::altExp()]) to switch for differential analysis.
#' The processing will be done in the main experiment while the differential
#' analysis will be done in the alternative experiment.
#' @param force_initial_clustering A logical specifying wether to force the 
#' initial number of cluster between 2 and 6. This is in order to avoid a too high
#' number of initial clusters which would be equivalent to a classical louvain
#' clustering.
#' @param verbose A logical specifying wether to print.
#'
#' @param ... 
#'
#' @return The SingleCellExperiment object with the assignation of cells to clusters.
#' If saving is true, also saves list of differential analyses, differential 
#' analyses summaries and embeddings for each re-clustered cluster. If runFDR is 
#' TRUE, also saves the list of FDR for each re-clusterd cluster.
#' 
#' 
#' @importFrom  Matrix Matrix rowSums
#' @importFrom qs qsave
#' @importFrom  SummarizedExperiment assays 
#' @importFrom  grDevices colors 
#' @importFrom  SingleCellExperiment reducedDim  
#' @import ggplot2
#' @import dplyr
#' @export 
#' @rdname iterative_differential_clustering
#' @exportS3Method iterative_differential_clustering default
#' 
iterative_differential_clustering.default <- function(
    object,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    n_dims = 10,
    dim_red = "PCA",
    vizualization_dim_red = "UMAP",
    processing_function = processing_ChromSCape,
    quantile.activation = 0.7,
    differential_function = differential_ChromSCape,
    logFC.th = log2(1.5),
    qval.th = 0.01,
    min_frac_cell_assigned = 0.1,
    limit = 5,
    starting.k = 100,
    starting.resolution = 0.1,
    resolution = 0.8,
    k = 50,
    FP_linear_model = NULL,
    color = NULL,
    nThreads = 10,
    swapExperiment = NULL,
    force_initial_clustering = TRUE,
    verbose = TRUE,
    ...
){
    set.seed(47)
    
    if(plotting == TRUE){
        if(dir.exists(output_dir)) 
            dir.create(file.path(output_dir, "iterations"), showWarnings = FALSE) else stop("output_dir is an invalid directory")
    }
    
    # For the first partition, try to find very low level clusters (e.g. low
    # resolution, high number of neighbors)
    if(force_initial_clustering){
      nclust = 1
      factor = 1
      max_iter = 10
      iter = 0
      while( (nclust < 2 | nclust > 6) & iter < max_iter){
        object. = ChromSCape::find_clusters_louvain_scExp(object,
                                                          k = starting.k,
                                                          resolution = factor * starting.resolution,
                                                          use.dimred = dim_red)
        nclust = length(unique(object.$cell_cluster))
        if(nclust < 2) factor = 1.15 * factor
        if(nclust > 6) factor = 0.85 * factor
        iter = iter + 1
      } 
      if(iter > max_iter) {
        warning("IDclust::iterative_differential_clustering - Didn't manage",
                " to find more than 1 initial cluster...")
        return()
      }
    } else {
      object. = ChromSCape::find_clusters_louvain_scExp(object,
                                                        k = starting.k,
                                                        resolution = starting.resolution,
                                                        use.dimred = dim_red)
    }

    object$IDcluster = gsub("C","A", object.$cell_cluster)
    rm(object.)
    gc()
    
    # Calculate the average % cells activated in a feature
    # and return a level of activation based on a given decile
    binmat = Matrix::Matrix((SummarizedExperiment::assays(object)$counts > 0) + 0, sparse = TRUE)
    min.pct.activation = quantile(Matrix::rowSums(binmat) / ncol(binmat), quantile.activation)
    
    # Find initial differences (even if there are none, the initial clusters
    # are always considered as true clusters).
    if(verbose) cat("Initial round of clustering - limit of differential genes set to 0",
                    " for this first round only.\n")
    
    if(length(swapExperiment)){
        if(verbose) cat("Swicthing experiment for DA to ", swapExperiment, ".\n")
      object = ChromSCape::swapAltExp_sameColData(object, alt = swapExperiment)
    }
    
    DA = find_differentiated_clusters(
        object, 
        differential_function = differential_function,
        by = "IDcluster",
        logFC.th = logFC.th,
        qval.th = qval.th,
        min_frac_cell_assigned = min_frac_cell_assigned,
        limit = 0,
        FP_linear_model = FP_linear_model,
        cluster_of_origin = "Omega",
        verbose = verbose,
        ...)
    
    if(length(swapExperiment)){
        if(verbose) cat("Swicthing back for clustering to main Experiment.\n")
        object = ChromSCape::getMainExperiment(object)
    }
    
    # Starting list of clusters to re-cluster
    differential_summary_df = DA$diffmat_n
    object$IDcluster = differential_summary_df$true_subcluster[match(object$IDcluster, differential_summary_df$subcluster)]
    
    
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
    object$cell_cluster_color = NULL
    
    # List of embeddings
    list_embeddings = list(SingleCellExperiment::reducedDim(object, dim_red))
    names(list_embeddings)[1] = paste(unique(object$IDcluster), collapse = "_")
    
    # List of marker features
    res = DA$res
    if(nrow(res) > 0) res$cluster_of_origin = "Omega"
    list_res = list(res)
    names(list_res)[1] = "Omega"
    
    iteration = 0
    gc()
    
    # Run IDC until no more clusters can be re-clustered into meaningful clusters
    while(iteration < nrow(differential_summary_df)){
        
        if(plotting == TRUE){
            # Plot each iteration of the algorithm
            png(file.path(output_dir, "iterations", paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
            object. = object
            object.$cell_cluster = gsub("Omega:","",object.$IDcluster)
            print(
                ChromSCape::plot_reduced_dim_scExp(object., reduced_dim = vizualization_dim_red, color_by = "IDcluster",
                                                   downsample = 50000, size = 0.35, transparency = 0.75, annotate_clusters = TRUE) +
                    ggtitle(paste0("Initital Clustering")) + theme(legend.position = "none")
            )
            dev.off()
            rm(object.)
        }
        
        iteration = iteration + 1
        
        if(verbose) cat("Doing iteration number ", iteration,"...\n")
        
        # Name of the cluster of origin
        partition_cluster_of_origin = differential_summary_df$true_subcluster[iteration]
        
        # Make sure that we do not re-cluster a 'parent cluster',
        # e.g. cells that were unassigned and therefore assigned to the parent
        # cluster
        if(!(partition_cluster_of_origin  %in% differential_summary_df$true_subcluster[duplicated(differential_summary_df$true_subcluster)])){
            
            # Letter assigned to the 'partition depth' (e.g. A, B, C...)
            partition_depth = which(LETTERS == substr(gsub(".*:","",partition_cluster_of_origin),1,1)) + 1
            
            # Select only cells from the given cluster
            object. = object[, which(object$IDcluster %in%  partition_cluster_of_origin)]
            
            if(verbose) cat("Re-calculating PCA and subclustering for cluster", partition_cluster_of_origin,
                            "with",ncol(object.),"cells.\n")
            
            # Re-cluster a cluster only if there are more than 100 cells
            if(ncol(object.) > 100){

                
                # Re-processing sub-cluster
                object. = processing_function(object., n_dims = n_dims, dim_red = dim_red)
                
                # Re-clustering sub-cluster
                object.. = ChromSCape::find_clusters_louvain_scExp(object., k = k, resolution =  resolution,
                                                                   use.dimred = dim_red)
                object.$IDcluster <- paste0(LETTERS[partition_depth],gsub("C", "", object..$cell_cluster))
                rm(object..)
                
                clusters = object.$IDcluster 
                cluster_u = unique(clusters)
                
                if(length(cluster_u) > 1 ){
                    if(verbose) cat("Found", length(cluster_u),"subclusters.\n")
                    if(verbose) print(table(clusters))
                    
                    
                    if(length(swapExperiment)){
                        if(verbose) cat("Swicthing experiment for DA to ", swapExperiment, ".\n")
                        object. = ChromSCape::swapAltExp_sameColData(object., alt = swapExperiment)
                    }
                    
                    # Find differentiated clusters
                    DA = find_differentiated_clusters(object.,
                                                      differential_function = differential_function,
                                                      by = "IDcluster",
                                                      logFC.th = logFC.th,
                                                      qval.th = qval.th,
                                                      min_frac_cell_assigned = min_frac_cell_assigned,
                                                      limit = limit,
                                                      FP_linear_model = FP_linear_model,
                                                      cluster_of_origin = partition_cluster_of_origin,
                                                      verbose = verbose,
                                                      ...)
                    gc()
                    
                    if(length(swapExperiment)){
                        if(verbose) cat("Swicthing back for clustering to main Experiment.\n")
                      object. = ChromSCape::getMainExperiment(object.)
                    }
                    
                    # Retrieve DA results
                    diffmat_n = DA$diffmat_n
                    res =  DA$res
                    res$cluster_of_origin = partition_cluster_of_origin[min(1,nrow(res))]
                    list_res[[partition_cluster_of_origin]] = res
                    list_embeddings[[partition_cluster_of_origin]] = SingleCellExperiment::reducedDim(object., dim_red)
                    
                    # If more than 'min_frac_cell_assigned' of the cells were assigned
                    # to 'true' subclusters (with marker features)
                    if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                        
                        # Add the new sublclusters to the list of clusters
                        differential_summary_df = rbind(differential_summary_df, diffmat_n)
                        object.$IDcluster = diffmat_n$true_subcluster[match(object.$IDcluster, diffmat_n$subcluster)]
                        object$IDcluster[match(colnames(object.), colnames(object))] = object.$IDcluster
                        
                        if(plotting == TRUE){
                            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                            print(
                                ChromSCape::plot_reduced_dim_scExp(object., color_by = "IDcluster",annotate_clusters = F,  reduced_dim = dim_red) + 
                                    ggtitle(paste0(partition_cluster_of_origin))
                            )
                            dev.off()
                        }
                    } else{
                        if(verbose) cat("Found 2 subcluster for", partition_cluster_of_origin," but with no difference. Maximum clustering reached...\n")
                    }
                } else{
                    if(verbose) cat("Found only 1 subcluster for", partition_cluster_of_origin,". Maximum clustering reached...\n")
                }
            }
            gc()
        } else{
            if(verbose) cat(partition_cluster_of_origin,"- This cluster was formed by 'unassigned' cells, not clustering it further...\n")
        }
    } 
    
    
    if(verbose) cat("\n\n\n##########################################################\nFinished !\nFound a total of", length(unique(object$IDcluster)),"clusters after",iteration ,"iterations.",
                    "\nThe average cluster size is ", floor(mean(table(object$IDcluster)))," and the median is",floor(median(table(object$IDcluster))),".",
                    "\nThe number of initital clusters not subclustered is ",length(grep(":", unique(object$IDcluster),invert = TRUE)),".",
                    "\n##########################################################\n")
    
    ## Saving results
    if(saving){
        # Table of differential features for each re-clustering
        IDC_DA = do.call("rbind", list_res)
        write.csv(IDC_DA, file = file.path(output_dir, "IDC_DA.csv"), quote = FALSE, row.names = FALSE)
        
        # List of embedding of each re-clustered cluster
        qs::qsave(list_embeddings, file.path(output_dir, "IDC_embeddings.qs"))
        
        # Summary table of the number of differential features for each re-clustering
        write.csv(differential_summary_df, file = file.path(output_dir, "IDC_summary.csv"), quote = FALSE, row.names = FALSE)
        
        # Final SingleCellExperiment with the clusters found by IDC 
        qs::qsave(object, file.path(output_dir, "scExp_IDC.qs"))
    }
    
    return(object)
}


#' @param object A Seurat object containing scRNA dataset with 'IDcluster' 
#' column.
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster. See [differential_edgeR_pseudobulk_LRT] for the default function.
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential 
#' passed to the differential_function.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential passed to the
#' differential_function.
#' @param limit An integer specifying the minimum number of features required 
#' for a subcluster to be called 'true' subcluster.
#' @param cluster_of_origin A character specifying the name of the cluster of 
#' origin that will be concatenated before the name of true subclusters. 
#' @param min_frac_cell_assigned A numeric between 0 and 1 specifying the
#' minimum percentage of the total cells in the SingleCellExperiment object that
#' needs to be assigned. If a lower proportion is assigned, all cells are 
#' assigned to the cluster of origin.
#' @param min_cluster_size An integer specifying the minimum number of cells
#' in a cluster to consider it as a 'true' subcluster.
#' @param verbose A logical specifying wether to print.
#' @param ... Additional parameters passed to the differential_function. See 
#' [differential_ChromSCape()] for more information on additional
#' parameters for the default function.
#' 
#' @return
#' @export 
#' @rdname find_differentiated_clusters
#' @exportS3Method find_differentiated_clusters Seurat
#' 
find_differentiated_clusters.Seurat <- function(object,
                                                differential_function = differential_edgeR_pseudobulk_LRT,
                                                by = "IDcluster",
                                                logFC.th = log2(1.5),
                                                qval.th = 0.01,
                                                limit = 5,
                                                cluster_of_origin = "Omega",
                                                min_frac_cell_assigned = 0.1,
                                                min_cluster_size = 30,
                                                verbose = TRUE,
                                                ...
){
    
    # If Seurat is a SingleCellExperiment or else, try convert it to Seurat
    if(!is(object, "Seurat")) {
        object = Seurat::as.Seurat(object)
        if(is.null(object$IDcluster) & !is.null(object$IDcluster)) 
            object$IDcluster = object$IDcluster else stop("Need 'IDcluster' or 'IDcluster', column.")
        Seurat::Idents(object) = object$IDcluster
    }
    
    set.seed(47)
    
    # Starting list of clusters to re-cluster
    cluster_u = as.character(unique(object@meta.data[[by]]))
    diffmat_n = data.frame(n_differential = rep(0,length(cluster_u)),
                           cluster_of_origin = cluster_of_origin,
                           subcluster = cluster_u,
                           true_subcluster = cluster_u)
    
    res = differential_function(object,
                                by = by,
                                logFC.th = logFC.th,
                                qval.th = qval.th,
                                ...)
    
    n_cell_assigned = 0
    for(i in seq_along(cluster_u)){
        group_cells = colnames(object)[which(object@meta.data[[by]] %in% cluster_u[i])]
        
        if(length(group_cells)  >= min_cluster_size){
            
            res. = res[which(res$cluster == cluster_u[i]),]
            diffmat_n$n_differential[i] = nrow(res.)
            if(verbose) cat(cluster_u[i],"- Found",diffmat_n$n_differential[i], "differential regions.\n")
            if(diffmat_n$n_differential[i] < limit){
                if(verbose) cat(cluster_u[i], " cluster has less than", limit, "enriched features.\nAssigning the cells to cluster of origin.\n")
                diffmat_n$true_subcluster[i] = cluster_of_origin
            } else{
                n_cell_assigned = n_cell_assigned + length(group_cells)
                diffmat_n$true_subcluster[i] = paste0(cluster_of_origin, ":", cluster_u[i])
            }
        } else {
            
        }
    }
    passing = TRUE
    
    if(verbose) cat("Finished finding differences - ", n_cell_assigned/ ncol(object), " fraction of cells were assigned.\n")
    if((n_cell_assigned/ ncol(object)) < min_frac_cell_assigned){
        if(verbose) cat("Not enough cells were assigned - not clustering.\n")
        passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}



#' @param object A Seurat object preprocessed with Seurat.
#' @param output_dir The output directory in which to plot and save objects.
#' @param plotting A logical specifying wether to save the plots or not.
#' @param saving A logical specifying wether to save the data or not.
#' @param limit An integer specifying the minimum number of significantly 
#' enriched / depleted features required in order for a subcluster to be called
#' a 'true' subcluster
#' @param dim_red The name of the slot to save the dimensionality reduction at
#' each step in the  \code{Seurat::Reductions(object)}.
#' @param vizualization_dim_red The name of the slot used for plotting. Must be
#' a valid slot present in \code{Seurat::Reductions(object)}.
#' @param processing_function A function that re-process the subset of clusters
#' at each step. It msut take in entry a Seurat object, \code{dim_red}
#'  and \code{n_dims} parameters and returns a Seurat containing a cell
#' embedding. See [processing_Seurat] for the default function.
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster. See [differential_edgeR_pseudobulk_LRT] for the default function.
#' @param logFC.th  A numeric specifying the log2 fold change of activation 
#' above/below which a feature is considered as significantly differential 
#' passed to the differential_function.
#' @param qval.th A numeric specifying the adjusted p-value below
#' which a feature is considered as significantly differential passed to the
#' differential_function.
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
#' @param n_dims An integer specifying the number of first dimensions to keep 
#' in the dimensionality reduction step.
#' @param nThreads  An integer specifying of threads to use
#' for the calculation of the FDR.
#' @param color Set of colors to use for the coloring of the clusters. This must
#' contains enough colors for each cluster (minimum 20 colors, but 100 colors
#' at least is recommended, based on the dataset).
#' @param force_initial_clustering A logical specifying wether to force the 
#' initial number of cluster between 2 and 6. This is in order to avoid a too high
#' number of initial clusters which would be equivalent to a classical louvain
#' clustering.
#' @param verbose A logical specifying wether to print.
#' @param ... Additional parameters passed to the differential_function. See 
#' [differential_edgeR_pseudobulk_LRT()] for more information on additional
#' parameters for the default function.
#'
#' @return The Seurat object with the assignation of cells to clusters.
#' If saving is true, also saves list of differential analyses, differential 
#' analyses summaries and embeddings for each re-clustered cluster.
#'
#' 
#' @import ggplot2
#' @import dplyr
#' @importFrom grDevices colors
#' @importFrom qs qsave
#' @export 
#' @rdname iterative_differential_clustering
#' @exportS3Method iterative_differential_clustering Seurat
#' 
iterative_differential_clustering.Seurat <- function(
    object,
    output_dir = "./",
    plotting = TRUE,
    saving = TRUE,
    n_dims = 50,
    dim_red = "pca",
    vizualization_dim_red = "umap",
    processing_function = processing_Seurat,
    differential_function = differential_edgeR_pseudobulk_LRT,
    logFC.th = log2(1.5),
    qval.th = 0.01,
    min_frac_cell_assigned = 0.1,
    limit = 10,
    starting.resolution = 0.1,
    starting.k = 100,
    resolution = 0.8,
    k = 100,
    color = NULL,
    nThreads = 10,
    force_initial_clustering = TRUE,
    verbose = TRUE,
    ...
){
    
    set.seed(47)
    
    if(plotting == TRUE){
        if(dir.exists(output_dir)) 
            dir.create(file.path(output_dir, "iterations"), showWarnings = FALSE) else stop("output_dir is an invalid directory")
    }
    
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
    object = Seurat::FindNeighbors(object, reduction = dim_red,  k.param = starting.k, dims = 1:n_dims, verbose = FALSE)
    
    if(force_initial_clustering){
      nclust = 1
      factor = 1
      max_iter = 10
      iter = 0
      while( (nclust < 2 | nclust > 6) & iter < max_iter){
        object = Seurat::FindClusters(object, algorithm = 2,
                                      resolution = factor * starting.resolution,
                                      random.seed = 47, verbose = FALSE)
        
        nclust = length(unique(object$seurat_clusters))
        if(nclust < 2) factor = 1.15 * factor
        if(nclust > 6) factor = 0.85 * factor
        iter = iter + 1
      } 
      if(iter > max_iter) {
        warning("IDclust::iterative_differential_clustering - Didn't manage",
                " to find more than 1 initial cluster...")
        return()
      }
    } else {
      object = Seurat::FindClusters(object, algorithm = 2, resolution = starting.resolution, random.seed = 47, verbose = FALSE)
      
    }
    object$IDcluster = object$seurat_clusters
    object$IDcluster <- paste0("A",as.numeric(object$IDcluster))
    Seurat::Idents(object) = object$IDcluster
    
    # Starting list of clusters to re-cluster
    cluster_u = unique(object$IDcluster)
    
    # Find initial differences (even if there are none, the initial clusters
    # are always considered as true clusters).
    if(verbose) cat("Initial round of clustering - limit of differential genes set to 0",
                    " for this first round only.\n")
    DA = find_differentiated_clusters(
        object,
        differential_function = differential_function,
        by = "IDcluster",
        logFC.th = logFC.th,
        qval.th = qval.th,
        min_frac_cell_assigned = min_frac_cell_assigned,
        cluster_of_origin =  "Omega",
        limit = 0,
        verbose = verbose,
        ...)
    
    # Starting list of clusters to re-cluster
    differential_summary_df = DA$diffmat_n
    object$IDcluster = differential_summary_df$true_subcluster[match(object$IDcluster, differential_summary_df$subcluster)]
    
    # List of embeddings
    list_embeddings = list(object@reductions[[dim_red]]@cell.embeddings )
    names(list_embeddings)[1] = paste(unique(object$IDcluster), collapse = "_")
    
    # List of differential analyses
    res = DA$res
    rownames(res) = NULL
    if(nrow(res) > 0)  res$cluster_of_origin = "Omega"
    list_res = list(res)
    names(list_res)[1] = "Omega"
    
    iteration = 0
    gc()
    
    while(iteration < nrow(differential_summary_df)){
        
        if(plotting == TRUE){
            # Plot initial
            png(file.path(output_dir, "iterations", paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
            object. = object
            object.$IDcluster = gsub("Omega:","",object.$IDcluster)
            print(
                Seurat::DimPlot(object, group.by =  "IDcluster",  reduction = vizualization_dim_red, cols = color) 
            )
            dev.off()
        }
        iteration = iteration + 1
        if(verbose) cat("Doing iteration number ", iteration,"...\n")
        
        
        partition_cluster_of_origin = differential_summary_df$true_subcluster[iteration]
        if(!(partition_cluster_of_origin  %in% differential_summary_df$true_subcluster[duplicated(differential_summary_df$true_subcluster)])){
            
            partition_depth = which(LETTERS == substr(gsub(".*:","",partition_cluster_of_origin),1,1)) + 1
            
            # Select first cluster
            object. = object[, which(object$IDcluster %in%  partition_cluster_of_origin)]
            if(verbose) cat("Re-calculating PCA and subclustering for cluster", partition_cluster_of_origin,".\n")
            
            if(ncol(object.) > 100){
                
                # Re-processing sub-cluster
                object. = processing_function(object., n_dims = n_dims, dim_red = dim_red)
                
                # Re-clustering sub-cluster
                object. = Seurat::FindNeighbors(object., reduction = dim_red, k.param = k, verbose = FALSE)
                object. = Seurat::FindClusters(object., algorithm = 2,   resolution = resolution,
                                               random.seed = 47, verbose = FALSE)
                object.$IDcluster = paste0(LETTERS[partition_depth], as.numeric(object.$seurat_clusters))
                Seurat::Idents(object.) = object.$IDcluster
                
                
                clusters = object.$IDcluster 
                cluster_u = unique(clusters)
                
                
                if(length(cluster_u) > 1 ){
                    if(verbose) cat("Found", length(cluster_u),"subclusters.\n")
                    if(verbose) print(table(clusters))
                    
                    ## Differential analysis
                    DA = find_differentiated_clusters(object., 
                                                      differential_function = differential_function,
                                                      by = "IDcluster",
                                                      logFC.th = logFC.th,
                                                      qval.th = qval.th,
                                                      min_frac_cell_assigned = min_frac_cell_assigned, 
                                                      limit = limit,
                                                      cluster_of_origin =  partition_cluster_of_origin,
                                                      verbose = verbose,
                                                      ...)
                    
                    # Retrieve DA results
                    diffmat_n = DA$diffmat_n
                    res = DA$res
                    rownames(res) = NULL
                    res$cluster_of_origin = partition_cluster_of_origin[min(1,nrow(res))]
                    list_res[[partition_cluster_of_origin]] = res
                    list_embeddings[[partition_cluster_of_origin]] = object.@reductions[[dim_red]]@cell.embeddings
                    
                    # If more than 'min_frac_cell_assigned' of the cells were assigned
                    # to 'true' subclusters (with marker features)
                    if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                        
                        # Add the new sublclusters to the list of clusters
                        differential_summary_df = rbind(differential_summary_df, diffmat_n)
                        object.$IDcluster = diffmat_n$true_subcluster[match(object.$IDcluster, diffmat_n$subcluster)]
                        object$IDcluster[match(colnames(object.), colnames(object))] = object.$IDcluster
                        
                        Seurat::Idents(object.) = object.$IDcluster
                        
                        if(plotting == TRUE){
                            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                            print(
                                Seurat::DimPlot(object., reduction = dim_red) + ggtitle(paste0(partition_cluster_of_origin, " - true"))
                            )
                            dev.off()
                        }
                    } else{
                        if(verbose) cat("Found 2 subcluster for",
                                        partition_cluster_of_origin," but with no difference. Maximum clustering reached...\n")
                    }
                    
                } else{
                    if(verbose) cat("Found only 1 subcluster for",
                                    partition_cluster_of_origin,". Maximum clustering reached...\n")
                }
            }
            gc()
            
        } else{
            if(verbose) cat(partition_cluster_of_origin,
                            "- This cluster was formed by 'unassigned' cells, not clustering it further...\n")
        }
    }
    
    if(verbose) cat("\n\n\n##########################################################\nFinished !\nFound a total of", length(unique(object$IDcluster)),"clusters after",iteration ,"iterations.",
                    "\nThe average cluster size is ",floor(mean(table(object$IDcluster)))," and the median is",floor(median(table(object$IDcluster))),".",
                    "\nThe number of initital clusters not subclustered is ",length(grep(":", unique(object$IDcluster),invert = TRUE)),".",
                    "\n##########################################################\n")
    
    
    ## Saving results
    if(saving){
        # Table of differential features for each re-clustering
        IDC_DA = do.call("rbind", list_res)
        write.csv(IDC_DA, file = file.path(output_dir, "IDC_DA.csv"), quote = FALSE, row.names = FALSE)
        
        # List of embedding of each re-clustered cluster
        qs::qsave(list_embeddings, file.path(output_dir, "IDC_embeddings.qs"))
        
        # Summary table of the number of differential features for each re-clustering
        write.csv(differential_summary_df, file = file.path(output_dir, "IDC_summary.csv"), quote = FALSE, row.names = FALSE)
        
        # Final SingleCellExperiment with the clusters found by IDC 
        qs::qsave(object, file.path(output_dir, "Seu_IDC.qs"))
    }
    
    return(object)
}

