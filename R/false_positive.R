#' Calculate FDR for SingleCellExperiment object
#'
#' @param object A SingleCellExperiment object 
#' @param by A character specifying the name of the metadata column referencing
#' the clusters.
#' @param differential_function A function that take in entry a 
#' SingleCellExperiment object and  parameters passed in ... and returns a 
#' data.frame containing the significantly differential features for each 
#' cluster.
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
calculate_FDR_scEpigenomics <- function(object,
                                        differential_function = differential_ChromSCape,
                                        by = "IDcluster",
                                        logFC.th = log2(1.5),
                                        qval.th = 0.01,
                                        nThreads = 1,
                                        iterations = 5,
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

        object$fake_cluster = sample(as.character(unlist(object[[by]])), length(as.character(unlist(object[[by]]))), replace = FALSE)
        df <- find_differentiated_clusters(
            object = object,
            differential_function = differential_function,
            by = "fake_cluster",
            logFC.th = logFC.th,
            qval.th = qval.th,
            min_frac_cell_assigned = 0,
            limit = 0,
            min_cluster_size = 0,
            verbose = FALSE,
            ...
            )
        df = df$diffmat_n
        df = df[,-c(2,4)]
        colnames(df) = c("n_differential", "cluster")
        df$iteration = i
        df
    }  
    
    FDR = do.call("rbind", FDR_iterations)
    
    parallel::stopCluster(myCluster)
    gc()
    
    return(FDR)
}