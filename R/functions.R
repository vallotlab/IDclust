geco.hclustAnnotHeatmapPlot.withColumn <- function (
    x = NULL, hc = NULL, hc_row = NULL, hmColors = NULL, 
    anocol = NULL, anorow = NULL, xpos = c(0.375, 0.9, 0.3745, 0.885, 0.05, 0.25),
    ypos = c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
    dendro.cex = 1, xlab.cex = 0.8, 
    hmRowNames = FALSE, hmRowNames.cex = 0.5, 
    hmColNames = FALSE, hmColNames.cex = 0.5,
    hmCategNamesRows = FALSE, hmCategNamesRows.cex = 0.5) 
{
    par(fig = c(xpos[1], xpos[2], ypos[5], ypos[6]), new = FALSE, 
        mar = c(0, 0, 1.5, 0))
    plot(hc, main = "", sub = "", las = 2, cex = dendro.cex, cex.axis = dendro.cex)
    par(fig = c(xpos[3], xpos[4], ypos[3], ypos[4]), new = TRUE, 
        mar = rep(0, 4))
    geco.imageCol(anocol, xlab.cex = xlab.cex, ylab.cex = 0)
    if (hmColNames) {
        axis(side = 3, padj = 0.5, hadj = 0.5,
             lwd = 0, at = seq(0, 1, length.out = ncol(x)), 
             labels = colnames(x), las = 2, cex.axis = hmColNames.cex)
    }
    par(fig = c(xpos[3], xpos[4], ypos[1], ypos[2]), new = TRUE, 
        mar = rep(0, 4))
    image(t(x), axes = FALSE, xlab = "", ylab = "", col = hmColors)
    par(fig = c(xpos[5], xpos[6], ypos[1] - 0.015, ypos[3] + 
                    0.015), new = TRUE, mar = c(0, 0, 0, 0))
    plot(as.phylo(hc_row), main = "", sub = "", las = 2, cex = dendro.cex, cex.axis = dendro.cex)
    par(fig = c(xpos[6] + 0.015, xpos[1], ypos[1], ypos[3]), 
        new = TRUE, mar = rep(0, 4))
    geco.imageCol(t(anorow)[, dim(anorow)[1]:1], xlab.cex = 0, 
                  ylab.cex = 0)
    box()
    if (hmRowNames) {
        axis(2, hadj = 0.1, lwd = 0, at = seq(0, 1, length.out = nrow(x)), 
             labels = rownames(x), las = 1, cex.axis = hmRowNames.cex)
    }
    if (hmCategNamesRows) {
        axis(1, hadj = 0.75, lwd = 0, at = seq(0, 1, length.out = ncol(anorow)), 
             labels = colnames(anorow), las = 2, cex.axis = hmCategNamesRows.cex)
    }
    
}

enrich_for_TF_ChEA3 <- function(genes_of_interest, n_random = 100, all_genes = unique(toTable(org.Hs.eg.db::org.Hs.egSYMBOL)[,2])){
    url = "https://maayanlab.cloud/chea3/api/enrich/"
    encode = "json"
    
    list_TF_enrichment = list()
    for(i in seq_len(n_random+1)){
        if(i != 1){
            genes = sample(all_genes, length(genes_of_interest))
        } else{
            genes = genes_of_interest
        }
        
        #POST to ChEA3 server
       
        httr::set_config(config(ssl_verifypeer = 0L))
        url = "https://maayanlab.cloud/chea3/api/enrich/"
        encode = "json"
        payload = list(query_name = "myQuery", gene_set = genes)
        
        #POST to ChEA3 server
        response = POST(url = url, body = payload, encode = encode)
        json = content(response, "text")
        
        #results as list of R dataframes
        results = fromJSON(json)
        
        #results as list of R dataframes
        list_TF_enrichment[[i]] = results$`Integrated--meanRank`
        
        if(i %% 10 == 0){cat("Done - ", i ," / ", n_random, ".\n")}
    }
    names(list_TF_enrichment) = c("Genes_of_interest", paste0("random_genes_",seq_len(n_random)))
    return(list_TF_enrichment)
}


TF.IDF.custom <- function(data, scale = 10000, log = TRUE, verbose = TRUE) {
    if (class(x = data) == "data.frame") {
        data <- as.matrix(x = data)
    }
    if (class(x = data) != "dgCMatrix") {
        data <- as(object = data, Class = "dgCMatrix")
    }
    if (verbose) {
        message("Performing TF-IDF normalization")
    }
    npeaks <- Matrix::colSums(x = data)
    tf <- t(x = t(x = data) / npeaks)
    # log transformation
    idf <- 1+ ncol(x = data) / Matrix::rowSums(x = data)
    norm.data <- Matrix::Diagonal(n = length(x = idf), x = idf) %*% tf
    if(log) norm.data = log1p(norm.data * scale) else norm.data = 
        norm.data * scale
    return(norm.data)
}


dice. = function (x, y) {
    M.11 = sum(x == 1 & y == 1)
    M.10 = sum(x == 1 & y == 0)
    M.01 = sum(x == 0 & y == 1)
    return (2*M.11 / (2*M.11 + M.10 + M.01))
}

dice <- function(m){
    dice_mat = matrix(0,ncol = nrow(m), nrow= nrow(m),
                      dimnames = list(rownames(m),
                                      rownames(m))) 
    for(row in 1:nrow(m)){
        for(col in 1:nrow(m)){
            dice_mat[row,col] = dice.(m[row,],m[col,])
        }
    }
    return(dice_mat)
}

jaccard. = function (x, y) {
    M.11 = sum(x == 1 & y == 1)
    M.10 = sum(x == 1 & y == 0)
    M.01 = sum(x == 0 & y == 1)
    return (M.11 / (M.11 + M.10 + M.01))
}

jaccard <- function(m){
    jac_mat = matrix(0,ncol = nrow(m), nrow= nrow(m),
                     dimnames = list(rownames(m),
                                     rownames(m))) 
    for(row in 1:nrow(m)){
        for(col in 1:nrow(m)){
            jac_mat[row,col] = jaccard.(m[row,],m[col,])
        }
    }
    return(jac_mat)
}

run_real_one_vs_all_comparisons_activation_scExp <- function(scExp,
                                                             name = "scExp",
                                                             qval.th = 0.1,
                                                             min.pct = 0.1,
                                                             logFC.th = 0.5,
                                                             limit = 5,
                                                             limit_by_proportion = NULL,
                                                             limit_factor = 2,
                                                             min.pct.cell_assigned = 0.75,
                                                             cluster_of_origin = "0"){
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
                    res[,paste0("group_activation.", group)] > min.pct
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
                             limit_by_proportion$mean_n_differential[index][1] + limit_factor * limit_by_proportion$sd_n_differential[index][1])
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
    if((n_cell_assigned/ ncol(scExp)) < min.pct.cell_assigned){
        cat("Not enough cells were assigned - not clustering.\n")
        passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}

run_real_pairwise_comparisons_Seurat <- function(Seu,
                                                 logFC.th = log2(1.5),
                                                 qval.th = 0.1,
                                                 test.use = c("Seurat", "edgeR-LRT"),
                                                 biological_replicate_col = "Target",
                                                 min.pct = 0.1,
                                                 limit = 10,
                                                 cluster_of_origin = "0",
                                                 min.pct.cell_assigned = 0.75){
  if(!is(Seu, "Seurat")) {
    Seu = as.Seurat(Seu)
    Seu$seurat_clusters = Seu$cell_cluster
    Idents(Seu) = Seu$seurat_clusters
    
  }
    cluster_u = unique(Seu$seurat_clusters)
    diffmat_n = data.frame(n_differential = rep(0,length(cluster_u)),
                           cluster_of_origin = cluster_of_origin,
                           subcluster = cluster_u,
                           subcluster_after_diff = cluster_u)
    if(test.use == "Seurat") {
      res = Seurat::FindAllMarkers(Seu, logfc.threshold = logFC.th, return.thresh = qval.th,
                                 test.use = "wilcox", min.pct = min.pct, only.pos = TRUE)
    } else {
      res = pseudobulk_DA_edgeR_Seurat(Seu, biological_replicate_col)
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
    if((n_cell_assigned/ ncol(Seu)) < min.pct.cell_assigned){
      cat("Not enough cells were assigned - not clustering.\n")
      passing = FALSE
    }
    
    out = list("diffmat_n" = diffmat_n, "res" = res, "passing_min_pct_cell_assigned" = passing)
    return(out)
}

pseudobulk_DA_edgeR_Seurat <- function(Seu, biological_replicate_col){
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
    y <- DGEList(counts=mat, group=group)
    keep <- filterByExpr(y)
    if(length(which(keep)) > 0){
      y <- y[keep,,keep.lib.sizes=FALSE]
      y <- calcNormFactors(y)
      design <- model.matrix(~group)
      y <- estimateDisp(y,design)
      fit <- glmFit(y,design)
      lrt <- glmLRT(fit,coef=2)
      topTags(lrt)
      tab = lrt$table
      
      binmat = Matrix((Seu@assays$RNA@counts > 0) + 0, sparse = TRUE)
      pct.1 = rowSums(binmat[,which(Seu$seurat_clusters == cluster_u[i])]) / length(which(Seu$seurat_clusters == cluster_u[i]))
      pct.2 = rowSums(binmat[,which(Seu$seurat_clusters != cluster_u[i])]) / length(which(Seu$seurat_clusters != cluster_u[i]))
      
      tab = tab %>% filter(abs(logFC) > 0.1 & PValue < 0.1) # very loose filter
      
      res. = data.frame(
        "p_val" = tab$PValue,
        "avg_log2FC"= tab$logFC,
        "pct.1" = pct.1[match(rownames(tab), names(pct.1))], 
        "pct.2"= pct.2[match(rownames(tab), names(pct.2))],
        "p_val_adj" = p.adjust(tab$PValue, method = "bonferroni"),
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

create_pseudobulk_mat_Seu <- function(Seu, biological_replicate_col){
  cluster_u = unique(Seu$seurat_clusters)
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
  mat = mat[,which(colSums(mat) > 0)]
  colnames(mat) = names_mat
  return(mat)
}
run_real_pairwise_comparisons_scExp_findMarkers <- function(scExp, logFC.th = log2(1.5), qval.th = 0.1,
                                                 test.use = "wilcox", min.pct = 0.1, limit = 10, name = "Seurat"){
    
    Seu = as.Seurat(scExp, counts = "counts", data = "normcounts")
    Seu$seurat_clusters = Seu$cell_cluster
    Idents(Seu) = Seu$seurat_clusters
    out = run_real_pairwise_comparisons_Seurat(Seu,  logFC.th = logFC.th, qval.th = qval.th,
                                               test.use = test.use, min.pct = min.pct, limit = limit, name = name)
    return(out)
}
run_real_pairwise_comparisons_scExp <- function(scExp, name = "scExp", logFC.th = log2(1.5), qval.th = 0.1, method = "wilcox", limit = 10, 
                                                BPPARAM = BiocParallel::bpparam(), plot = FALSE, output = "./"){
    
    cluster_u = unique(names(sort(table(scExp$cell_cluster))))
    diffmat_n = matrix(0, nrow = length(cluster_u), ncol = length(cluster_u), dimnames = list(cluster_u,cluster_u))
    res_list = list()
    for(i in 1:(length(cluster_u) -1) ){
        for(j in (i+1):length(cluster_u) ){
            cat("Differential analysis", cluster_u[i],"vs", cluster_u[j],"...")
            group_df = setNames(data.frame(cluster_u[i]), cluster_u[i])
            ref_df = setNames(data.frame(cluster_u[j]), cluster_u[j])
            scExp. = scExp[,which(scExp$cell_cluster %in% c(cluster_u[i],cluster_u[j]))]
            scExp. = scExp.[which(rowSums(counts(scExp.)) > quantile(rowSums(counts(scExp.)),0.1)),]
            scExp. = differential_analysis_scExp(scExp., de_type = "custom", method = method,
                                                 qval.th = qval.th, cdiff.th = logFC.th, prioritize_genes = F,
                                                 group = group_df, ref = ref_df, BPPARAM =BPPARAM)
            res = scExp.@metadata$diff$res
            diffmat_n[i,j] = diffmat_n[j,i] = scExp.@metadata$diff$summary[1]
            
            cat("found", scExp.@metadata$diff$summary[1], "differential regions.\n")
            if(scExp.@metadata$diff$summary[1] < limit){
                cat("Exiting because one cluster had less than", limit, "differential regions.\n")
                return(FALSE)
            }
            if(plot){
                png(file.path(output, paste0(name,"_volcano_",cluster_u[i],"_vs_",cluster_u[j],".png")), width = 1200, height = 1200, res = 200)
                print(
                    ChromSCape::plot_differential_volcano_scExp(scExp, cell_cluster = cluster_u[i], cdiff.th = logFC.th, qval.th =qval.th)
                )
                dev.off()
             }
            
            cluster_name = paste0(name,":",cluster_u[i],"_vs_", cluster_u[j]) 
            res_list[[cluster_name]] = scExp.@metadata$diff$res
        }
    }
    out = list("diffmat_n" = diffmat_n, "res_list" = res_list)
    return(out)
}
# 
# scExp = qs::qread("output/Marks/PairedTag/H3K4me1/scExp_5000.qs")
# cluster = "A1"
# cluster_mat = counts(scExp[,which(scExp$cell_cluster == "A1")])
# reference_mat = counts(scExp[,which(scExp$cell_cluster == "A2")])
# reference_mat = Matrix((reference_mat> 0) + 0, sparse = TRUE)
# is_binned_reference_mat = TRUE
# plot = TRUE
# name = "cluster"
# differential_activation <- function(cluster_mat, reference_mat, is_binned_reference_mat = FALSE,
#                                     plot = TRUE, name = "cluster"){
#     cluster_bin_mat =  Matrix((cluster_mat> 0) + 0, sparse = TRUE)
#     rectifier = mean(colSums(cluster_mat)) / mean(colSums(reference_mat))
#     group_sum = Matrix::rowSums(cluster_bin_mat)
#     group_activation = group_sum / (ncol(cluster_bin_mat) * rectifier)
# 
#     reference_bin_mat = reference_mat
#     if(!is_binned_reference_mat) reference_bin_mat = Matrix((reference_bin_mat> 0) + 0, sparse = TRUE)
# 
#     reference_sum = Matrix::rowSums(reference_bin_mat) 
#     reference_activation = reference_sum / ncol(reference_bin_mat)
# 
# 
#     n_cell_cluster = ncol(cluster_bin_mat)
#     n_cell_reference = ncol(reference_bin_mat)
#     other_group =  n_cell_cluster- group_sum
#     other_ref = n_cell_reference - reference_sum
#     
#     new_env = new.env()
#     new_env$chisq_mat = cbind(group_sum, other_group, reference_sum, other_ref)
#     
#     cl <- parallel::makeCluster(getOption("cl.cores", 6))
#     parallel::clusterExport(cl, varlist = c("chisq_mat"), envir = new_env)
#     pvalues = parallel::parApply(cl = cl, new_env$chisq_mat, 1, 
#                                  function(row) chisq.test(matrix(c(row[1], row[2], row[3], row[4]), ncol = 2),simulate.p.value = FALSE)$p.value)
#     rm(cl)
#     
#     q.values = p.adjust(pvalues, method = "BH")
#     logFCs = log2(group_activation/reference_activation)
#    
#     df = data.frame(features = rownames(cluster_mat),
#                     group_activation = group_activation,
#                     real_group_percentage = group_sum / ncol(cluster_bin_mat),
#                     reference_activation = reference_activation,
#                     logFC = logFCs,
#                     q.value = q.values)
#     df = df[order(df$q.value, df$logFC, decreasing = F),]
#     if(plot){
#         
#         h = hist(reference_activation, col = alpha("red",0.2), breaks = 75,
#                  main = paste0(name," vs the rest - rectifier = ", round(rectifier,2)), xlab = "% Activation in reference")
#         # par(new=T)
#         # plot(density(negbin_distribution / ncol(reference_bin_mat)), col = "red", axes=F,main=NA,xlab=NA,ylab=NA, xlim = c(0, max(h$breaks)))
#         par(new=T)
#         hist( group_sum / (ncol(cluster_bin_mat)), col = alpha("dark green",0.2),breaks = 75, xlim = c(0, max(h$breaks)), axes=F,main=NA,xlab=NA,ylab=NA)
#         
#         # top = head(df$features, min(nrow(df), 5))
#         # abline(v=reference_activation[top], col = "black", lty = 3)
#         # abline(v=group_activation[top], col = "blue", lty = 3)
#         # 
#         # quant_thresh = quantile(poisson_distribution/(ncol(reference_bin_mat)), 1 - qval.th)
#         # abline(v=quant_thresh, lty = 3, col = "red")
#         # text(x=quant_thresh+0.01, y = max(h$counts) / 2, labels = paste0("Q.",names(quant_thresh)), col = "red")
#         
#         legend((2/5) * range(h$breaks)[2], y = nrow(reference_bin_mat)/ 15,
#                legend=c("Reference Distribution", "Group Distribution"),
#                col=c("red", "dark green"), lty=1:3, cex=0.8)
#     }
#     
#     return(df)
# }

find_two_clusters_louvain_Seurat <- function(Seu, starting_nclust = length(scExp$cell_cluster),
                                                   starting_resolution = 0.8,  use.dimred = "PCA",
                                                   incrementation  = 0.5, nPCA = 50){
    n_clust = starting_nclust
    resolution = starting_resolution
    while(n_clust > 2){
        resolution = resolution * incrementation
        cat("Clustering with resolution = ", resolution, "\n")
        clusters = find_clusters_louvain_Seurat(Seu, resolution = resolution, use.dimred = use.dimred, nPCA = nPCA) 
        n_clust = length(unique(clusters))
        cat("Found", n_clust, "clusters with resolution =", resolution,".\n")
    }
    return(clusters)
}

find_two_clusters_louvain_Seurat_scExp <- function(scExp, starting_nclust = length(scExp$cell_cluster),
                                            starting_resolution = 0.8, k = 100,  use.dimred = "PCA",
                                            incrementation  = 0.5, nPCA = 50){
    n_clust = starting_nclust
    resolution = starting_resolution
    while(n_clust > 2){
        resolution = resolution * incrementation
        cat("Clustering with resolution = ", resolution, "\n")
        clusters = find_clusters_louvain_Seurat_scExp(scExp, k = k,  resolution = resolution, use.dimred = use.dimred, nPCA = nPCA) 
        n_clust = length(unique(clusters))
        cat("Found", n_clust, "clusters with k =", k,".\n")
    }
    return(clusters)
}

find_two_clusters_louvain_scExp <- function(scExp, starting_nclust = length(scExp$cell_cluster),
                                            starting_k =100, type ="jaccard", use.dimred = "PCA",
                                            incrementation  = 0.5 ){
    n_clust = starting_nclust
    k = starting_k
    while(n_clust > 2){
        k = (1 + incrementation) * k
        cat("Clustering with k = ", k, "\n")
        clusters = find_clusters_louvain_scExp(scExp, k = k, type = type, use.dimred = use.dimred) 
        n_clust = length(unique(clusters))
        cat("Found", n_clust, "clusters with k =", k,".\n")
    }
    return(clusters)
}

find_clusters_louvain_Seurat <- function(Seu,  resolution = 0.8, k = 20, use.dimred = "pca", nPCA = 10){
    Seu = Seurat::FindNeighbors(Seu, reduction = use.dimred,  k.param = k, dims = 1:nPCA, verbose =F)
    Seu = Seurat::FindClusters(Seu, algorithm = 2, resolution = resolution, random.seed = 47, verbose = FALSE)
    return(as.numeric(Seu$seurat_clusters))
}

find_clusters_louvain_Seurat_scExp <- function(scExp, resolution = 0.8,  k = 20, use.dimred = "PCA", nPCA = 10){
    seu = as.Seurat(scExp, counts = "counts", data = "normcounts")
    clusters = find_clusters_louvain_Seurat(seu, resolution = resolution, k = k, use.dimred = use.dimred, nPCA = nPCA)
    return(paste0("C", as.character(as.numeric(clusters))))
}


import_scExp_gz <- function(datadir, pattern){
    files = list.files(datadir, pattern = pattern, full.names = TRUE)
    dir.create(file.path(datadir,"tmp"))
    files = sapply(files, function(f){
        tmp = file.path(datadir,"tmp",paste0(gsub(pattern,"",basename(f)),".tsv"))
        gunzip(f, tmp, overwrite = TRUE, remove = FALSE)
    })
    out <- ChromSCape::import_scExp(files)
    unlink(file.path(datadir,"tmp"),recursive = TRUE,force =  TRUE)
    return(out)
}

plot_reduced_dim_scExp_devel <- function(scExp, color_by = "sample_id", reduced_dim = c("PCA", 
                                                                                        "TSNE", "UMAP"),
                                         select_x = "Component_1",
                                         select_y = "Component_2",
                                         downsample = 5000,
                                         transparency = 0.6,
                                         size = 1)
{
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.character(color_by), is.character(reduced_dim), 
              is.character(select_x), is.character(select_y), is.numeric(5000),
              is.numeric(transparency))
    
    if (!reduced_dim[1] %in% SingleCellExperiment::reducedDimNames(scExp)) 
        stop(paste0("ChromSCape::plot_reduced_dim_scExp - ", reduced_dim[1], " is not present in object, please run normalize_scExp first."))
    
    if (!color_by %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by must be present in colnames of colData(scExp).")
    
    if (!paste0(color_by, "_color") %in% colnames(SingleCellExperiment::colData(scExp))) 
        stop("ChromSCape::plot_reduced_dim_scExp - color_by's color column must be present 
         in colnames of colData(scExp). Please run color_scExp first.")
    
    if (!select_x %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_x must be present in colnames of PCA of scExp.")
    
    if (!select_y %in% colnames(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]))) 
        stop("ChromSCape::plot_reduced_dim_scExp - select_y must be present in colnames of PCA of scExp.")
    
    set.seed(2047)
    if(ncol(scExp) > downsample ) scExp = scExp[,sample(ncol(scExp),downsample,replace = F)]
    
    plot_df = as.data.frame(cbind(SingleCellExperiment::reducedDim(scExp, reduced_dim[1]), 
                                  SingleCellExperiment::colData(scExp)))
    
    p <- ggplot(plot_df, aes_string(x = select_x, y = select_y)) + geom_point(alpha = transparency,
                                                                              size = size, aes(color = SingleCellExperiment::colData(scExp)[, color_by])) +
        labs(color = color_by) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill = NA))
    
    if (color_by == "total_counts")
    {
        p <- p + scale_color_gradientn(colours = matlab.like(100))
    } else
    {
        
        cols = unique(as.character(
            SingleCellExperiment::colData(scExp)[,paste0(color_by, "_color")]))
        names(cols) = unique(as.character(
            SingleCellExperiment::colData(scExp)[,color_by]))
        p <- p + scale_color_manual(values = cols)
    }
    return(p)
}

makeAttr <- function(graph, default, valNodeList) {
    tmp <- nodes(graph)
    x <- rep(default, length(tmp)); names(x) <- tmp
    
    if(!missing(valNodeList)) {
        stopifnot(is.list(valNodeList))
        allnodes <- unlist(valNodeList)
        stopifnot(all(allnodes %in% tmp))
        for(i in seq(valNodeList)) {
            x[valNodeList[[i]]] <- names(valNodeList)[i]
        }
    }
    return(x)
}

show_in_excel <- function(.data){
    tmp <- paste0(tempfile(),".csv")
    write.csv(.data, tmp)
    fs::file_show(path = tmp)
}

VennDiagram_3 <- function(genesSet1, genesSet2,genesSet3,  savePNG=F,png_title,title,groups) {
    genesSet1  = as.vector(genesSet1)
    genesSet2  = as.vector(genesSet2)
    genesSet3  = as.vector(genesSet3)
    
    aera1 = length(genesSet1)
    aera2 = length(genesSet2)
    aera3 = length(genesSet3)
    
    n12 = length(intersect(genesSet1,genesSet2))
    n23 = length(intersect(genesSet2,genesSet3))
    n13 = length(intersect(genesSet1,genesSet3))
    n123 = length(intersect(intersect(genesSet1,genesSet2),genesSet3))
    
    if(savePNG == TRUE){
        png(paste(png_title,title,".png",sep=""))
        grid.newpage()
        draw.triple.venn(area1 = aera1, area2 = aera2, area3 = aera3, n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = groups, 
                         fill = c("light green", "red", "orange"),euler.d=T, scaled=T)
        dev.off()
    }
    else{
        grid.newpage()
        draw.triple.venn(area1 = aera1, area2 = aera2, area3 = aera3, n12 = n12, n23 = n23, n13 = n13, n123 = n123, category = groups, 
                         fill = c("light green", "red", "orange"),euler.d = T, scaled=T)
    }
}

VennDiagram_2 <- function(genesSet1, genesSet2, savePNG=F,png_title="VennDiagram",title="Main",groups) {
    genesSet1  = as.vector(genesSet1)
    genesSet2  = as.vector(genesSet2)
    aera1 = length(genesSet1)
    aera2 = length(genesSet2)
    
    
    cross.area = length(intersect(genesSet1,genesSet2))
    
    if(savePNG == TRUE){
        png(paste(png_title,title,".png",sep=""))
        
        grid.newpage()
        draw.pairwise.venn(area1 = aera1, area2 = aera2, cross.area = cross.area, category = groups,
                           , fill = c("light green", "red"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
        dev.off()
    }
    else{
        grid.newpage()
        draw.pairwise.venn(area1 = aera1, area2 = aera2, cross.area = cross.area, category =groups,
                           , fill = c("light green", "red"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
    }
}

violin_plot_cooccurence <- function(
    mat_init,mat_pers = NULL,
    sub_cells = list(colnames(mat_pers),colnames(mat_init)),
    levels=c("Untreated","Persister"),
    colors =c("#afafafff", "#1AB8AD")){
    
    if(!is.null(mat_pers)) mat = Matrix::cBind(mat_pers,mat_init)
    else mat = mat_init
    
    # mat = mat[which(rownames(mat) %in% overexpressed_MM468_persister_regions),]
    mat = mat[,which(colnames(mat) %in% unlist(sub_cells))]
    bin_mat = mat
    bin_mat[(bin_mat>0)]=1
    
    coocurrence_persister_genes_score = as.data.frame(Matrix::colSums(bin_mat) / nrow(bin_mat))
    colnames(coocurrence_persister_genes_score) = "coocurrence_score"
    rownames(coocurrence_persister_genes_score) = colnames(bin_mat)
    coocurrence_persister_genes_score$sample = ""
    coocurrence_persister_genes_score[sub_cells[[1]],"sample"] = levels[1]
    if(!is.null(mat_pers)) coocurrence_persister_genes_score[sub_cells[[2]],"sample"] = levels[2]
    
    coocurrence_persister_genes_score$sample = factor(
        coocurrence_persister_genes_score$sample,levels=levels)
    
    p = ggplot(coocurrence_persister_genes_score, aes(x=sample,y=coocurrence_score,
                                                      fill=sample)) + 
        geom_violin(alpha=0.8) + theme_classic() + scale_fill_manual(values=colors) +
        stat_summary(fun=median, geom="point", size=2, color="black") + geom_jitter(size=0.2, alpha=0.2) 
    if(!is.null(mat_pers)) p = p + ggpubr::stat_compare_means(method = "t.test",ref.group = levels[1])
    p
}
