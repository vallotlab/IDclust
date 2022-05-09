
scExp_differential_clustering <- function(
    scExp,
    starting_differential_summary_df,
    output_dir = "./",
    nPCA = 10,
    percent_feature = 1,
    quantile.activation = 0.75,
    min.pct.cell_assigned = 0.25,
    FC.th = 5,
    qval.th = 0.1,
    limit = 5,
    k = 100,
    resolution = 0.1,
    runFDR = TRUE,
    limit_by_proportion = NULL,
    nThreads = 10
){
    FDR_list = list()
    set.seed(47)
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
    color <- sample(color, 400)
    scExp$cell_cluster_color = NULL
    differential_summary_df = starting_differential_summary_df
    
    list_embeddings = list(reducedDim(scExp,"PCA"))
    names(list_embeddings)[1] = paste(unique(scExp$cell_cluster), collapse = "_")
    
    list_res = list()
    
    diffmat_list = list()
    iteration = 0
    gc()
    
    # calculate the average % cells activated in a feature
    # and return a level of activation based on a given decile
    binmat = Matrix((assays(scExp)$counts > 0) + 0, sparse = TRUE)
    min.pct.activation = quantile(Matrix::rowSums(binmat) / ncol(binmat), quantile.activation)
    
    
    while(iteration < nrow(differential_summary_df)){
        
        # Plot initial
        png(file.path(output_dir, paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
        print(
            plot_reduced_dim_scExp(scExp, reduced_dim = "UMAP", color_by = "cell_cluster", downsample = 50000, size = 0.35, transparency = 0.75) +
                ggtitle(paste0("Initital Clustering")) + theme(legend.position = "none")
        )
        dev.off()
        
        iteration = iteration + 1
        cat("Doing iteration number ", iteration,"...\n")
        
        partition_cluster_of_origin = differential_summary_df$cluster[iteration]
        if(!(partition_cluster_of_origin  %in% differential_summary_df$cluster[duplicated(differential_summary_df$cluster)])){
        partition_depth = which(LETTERS == substr(gsub(".*:","",partition_cluster_of_origin),1,1)) + 1
        
        # Select first cluster
        scExp. = scExp[, which(scExp$cell_cluster %in%  partition_cluster_of_origin)]
        
        cat("Re-calculating PCA and subclustering for cluster", partition_cluster_of_origin,
            "with",ncol(scExp.),"cells.\n")
        cat()
        if(ncol(scExp.) > 100){
            
            # Rerun PCA
            scExp. = ChromSCape::find_top_features(scExp.,n = floor(percent_feature * nrow(scExp.)), keep_others = FALSE)
            scExp. = ChromSCape::normalize_scExp(scExp., type = "TFIDF")
            scExp. = ChromSCape::reduce_dims_scExp(scExp., dimension_reductions = "PCA", n = nPCA, remove_PC = "Component_1")
            scExp.$cell_cluster = find_clusters_louvain_Seurat_scExp(scExp., resolution =  resolution, use.dimred = "PCA", nPCA = nPCA-1)
            scExp.$cell_cluster <- paste0(LETTERS[partition_depth],gsub("C", "", scExp.$cell_cluster))
            
            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_raw.png")), width = 1400, height = 1200, res = 200)
            print(
                plot_reduced_dim_scExp(scExp., color_by = "cell_cluster", reduced_dim = "PCA") + 
                    ggtitle(paste0(partition_cluster_of_origin, " - raw"))
            )
            dev.off()
            clusters = scExp.$cell_cluster 
            cluster_u = unique(clusters)
            
            if(length(cluster_u) > 1 ){
                cat("Found", length(cluster_u),"subclusters.\n")
                print(table(clusters))
                

                gc()
                ## Differential analysis
                DA = run_real_one_vs_all_comparisons_activation_scExp(scExp = scExp.,
                                                                      qval.th = qval.th,
                                                                      min.pct = min.pct.activation,
                                                                      logFC.th = log2(FC.th),
                                                                      min.pct.cell_assigned = min.pct.cell_assigned,
                                                                      limit = limit,
                                                                      limit_by_proportion = limit_by_proportion,
                                                                      cluster_of_origin = partition_cluster_of_origin)
                gc()
                diffmat_n = DA$diffmat_n
                print(diffmat_n)
                list_res[[partition_cluster_of_origin]] = DA$res
                list_embeddings[[partition_cluster_of_origin]] = reducedDim(scExp., "PCA")
                diffmat_list[[partition_cluster_of_origin]] = diffmat_n
                
                
               if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                        
                        df = data.frame("cluster" = diffmat_n$subcluster_after_diff, n_differential = diffmat_n$n_differential)
                        differential_summary_df = rbind(differential_summary_df, df)
                        scExp.$cell_cluster = diffmat_n$subcluster_after_diff[match(scExp.$cell_cluster, diffmat_n$subcluster)]
                        scExp$cell_cluster[match(colnames(scExp.), colnames(scExp))] = scExp.$cell_cluster

                    if(runFDR == TRUE){
                        gc()
                        myCluster <- makeCluster(nThreads, type = "FORK")
                        doParallel::registerDoParallel(myCluster)
                        scExp.. = scExp.
                        FDR_list. = foreach(i = 1:10) %dopar% {
                            set.seed(i * 47)
                            scExp..$cell_cluster = sample(scExp..$cell_cluster, length(scExp..$cell_cluster), replace = F)
                            DA <- run_real_one_vs_all_comparisons_activation_scExp(scExp = scExp..,
                                                                             qval.th = qval.th,
                                                                             min.pct = min.pct.activation,
                                                                             logFC.th = log2(FC.th),
                                                                             min.pct.cell_assigned = min.pct.cell_assigned,
                                                                             limit = limit,
                                                                             cluster_of_origin = partition_cluster_of_origin)
                            df = DA$diffmat_n[,-4]
                            df$iteration = i
                            df
                        }  
                        FDR_list[[partition_cluster_of_origin]] = do.call("rbind", FDR_list.)
                        FDR = length(which(FDR_list[[partition_cluster_of_origin]]$n_differential >= 5)) / nrow(FDR_list[[partition_cluster_of_origin]])
                        parallel::stopCluster(myCluster)
                        gc()
                    }
                        
                        tit = ifelse(runFDR, paste0(" - FDR - ", FDR), "")
                        png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                        print(
                            plot_reduced_dim_scExp(scExp., color_by = "cell_cluster", reduced_dim = "PCA") + 
                                ggtitle(paste0(partition_cluster_of_origin, tit))
                        )
                        dev.off()
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
    
    qs::qsave(list_res, file.path(output_dir, "list_res.qs"))
    qs::qsave(list_embeddings, file.path(output_dir, "list_embeddings.qs"))
    qs::qsave(diffmat_list, file.path(output_dir, "diffmat_list.qs"))
    qs::qsave(scExp, file.path(output_dir, "scExp_clustered.qs"))
    qs::qsave(FDR_list, file.path(output_dir, "FDR_list.qs"))
    
}

