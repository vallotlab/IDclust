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
                                        by = "IDcluster",
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
        scExp$IDcluster = sample(scExp$IDcluster, length(scExp$IDcluster), replace = FALSE)
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