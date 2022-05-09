source("scripts/functions.R")
Seurat_differential_clustering <- function(
    Seu,
    method = c("Seurat", "monocle")[1],
    starting_differential_summary_df,
    output_dir = "./",
    nPCA = 50,
    nfeatures = 2000,
    logFC.th = log2(2),
    qval.th = 0.01,
    min.pct = 0.1,
    min.pct.cell_assigned = 0.25,
    limit = 10,
    test.use = "wilcox",
    k = 100,
    resolution = 0.5,
    biological_replicate_col = "Target"
){
    
    set.seed(47)
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
    color <- sample(color, 400)

    differential_summary_df = starting_differential_summary_df
    
    list_embeddings = list(Seu@reductions$pca@cell.embeddings )
    names(list_embeddings)[1] = paste(unique(Seu$seurat_clusters), collapse = "_")
    
    list_res = list()
    
    diffmat_list = list()
    iteration = 0
    gc()

    while(iteration < nrow(differential_summary_df)){
        
        # Plot initial
        png(file.path(output_dir, paste0("Iteration_",iteration,".png")), width = 1600, height = 1200, res = 200)
        print(
            DimPlot(Seu, group.by =  "seurat_clusters", cols = color) 
        )
        dev.off()
        
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
            Seu. = FindVariableFeatures(Seu., selection.method = "vst", verbose = FALSE)
            Seu. <- ScaleData(Seu., verbose = FALSE)
            Seu. <- RunPCA(Seu., features = VariableFeatures(object = Seu.), verbose = FALSE)
            Seu. <- FindNeighbors(Seu., use.dimred = "pca", verbose = FALSE)
            Seu. <- FindClusters(Seu., algorithm = 2, resolution = 0.1,  random.seed = 47, verbose = FALSE)
            Seu.$seurat_clusters <- paste0(LETTERS[partition_depth],as.numeric(Seu.$seurat_clusters))
            Idents(Seu.) = Seu.$seurat_clusters
            
            
            png(file.path(output_dir, paste0(partition_cluster_of_origin,"_raw.png")), width = 1400, height = 1200, res = 200)
            print(
                DimPlot(Seu., reduction = "pca") + ggtitle(paste0(partition_cluster_of_origin, " - raw"))
            )
            dev.off()
            
            clusters = Seu.$seurat_clusters 
            cluster_u = unique(clusters)
            
            
            if(length(cluster_u) > 1 ){
                cat("Found", length(cluster_u),"subclusters.\n")
                print(table(clusters))
                
                ## Differential analysis
                DA = run_real_pairwise_comparisons_Seurat(Seu., logFC.th = logFC.th, qval.th = qval.th,
                                                          min.pct = min.pct, test.use = test.use, biological_replicate_col = biological_replicate_col,
                                                          min.pct.cell_assigned = min.pct.cell_assigned, 
                                                          cluster_of_origin =  partition_cluster_of_origin, limit = limit)
                
                diffmat_n = DA$diffmat_n
                print(diffmat_n)
                list_res[[partition_cluster_of_origin]] = DA$res
                list_embeddings[[partition_cluster_of_origin]] = Seu.@reductions$pca@cell.embeddings
                diffmat_list[[partition_cluster_of_origin]] = diffmat_n
                
                    if(!isFALSE(DA$passing_min_pct_cell_assigned)){
                      df = data.frame("cluster" = diffmat_n$subcluster_after_diff, n_differential = diffmat_n$n_differential)
                      differential_summary_df = rbind(differential_summary_df, df)
                      
                      Seu.$seurat_clusters = diffmat_n$subcluster_after_diff[match(Seu.$seurat_clusters, diffmat_n$subcluster)]
                      Seu$seurat_clusters[match(colnames(Seu.), colnames(Seu))] = Seu.$seurat_clusters
                      
                      list_embeddings[[partition_cluster_of_origin]] = Seu.@reductions$pca@cell.embeddings
                      diffmat_list[[partition_cluster_of_origin]] = diffmat_n
                      
                      Idents(Seu.) = Seu.$seurat_clusters
                      png(file.path(output_dir, paste0(partition_cluster_of_origin,"_true.png")), width = 1400, height = 1200, res = 200)
                      print(
                          DimPlot(Seu., reduction = "pca") + ggtitle(paste0(partition_cluster_of_origin, " - true"))
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
    
    cat("\n\n\n##########################################################\nFinished !\nFound a total of", length(unique(Seu$seurat_clusters)),"clusters after",iteration ,"iterations.",
        "\nThe average cluster size is ",floor(mean(table(Seu$seurat_clusters)))," and the median is",floor(median(table(Seu$seurat_clusters))),".",
        "\nThe number of initital clusters not subclustered is ",length(grep(":", unique(Seu$seurat_clusters),invert = TRUE)),".",
        "\n##########################################################\n")
    
    qs::qsave(list_res, file.path(output_dir, "list_res.qs"))
    qs::qsave(list_embeddings, file.path(output_dir, "list_embeddings.qs"))
    qs::qsave(diffmat_list, file.path(output_dir, "diffmat_list.qs"))
    qs::qsave(Seu, file.path(output_dir, "Seu_clustered.qs"))
    
}

