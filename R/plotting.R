#' Plot Iterative Differential Clustering network
#'
#' @param object  An object clustered with
#'  [iterative_differential_clustering()].
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' # Plotting of Seurat scRNA object (Paired-Tag)
#' if(requireNamespace("Seurat", quietly=TRUE)){
#' 
#' data("Seu", package = "IDclust")
#' data("IDC_summary_scRNA", package = "IDclust")
#' 
#' plot_cluster_network(
#'     object = Seu,
#'     IDC_summary = IDC_summary_scRNA,
#'     color_by = "IDcluster",
#'     cluster_col = "IDcluster",
#'     colors = NULL,
#'     node_size_factor = 7.5,
#'     edge_size_factor = 1,
#'     function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE,
#'                                                          flip.y = FALSE)
#' )
#' 
#' # Plotting proportion of cells activating a specific gene in  Seurat scRNA 
#' # object (Paired-Tag)
#' plot_cluster_network(
#'     object = Seu,
#'     IDC_summary = IDC_summary_scRNA,
#'     color_by = "Erbb4", # a gene contained in the Seu object
#'     threshold_to_define_feature_active = 2,
#'     assay = "RNA",
#'     cluster_col = "IDcluster",
#'     colors = NULL,
#'     node_size_factor = 7.5,
#'     edge_size_factor = 1,
#'     function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE,
#'                                                          flip.y = FALSE)
#' )
#' 
#' # Plotting cluster networ with the pathway information on the edges:
#' data("IDC_DA_scRNA", package = "IDclust")
#' edge_df = top_enriched_pathways(
#'     IDC_DA_scRNA,
#'     top = 5,
#'     gene_col = "gene",
#'     qval.th = 0.1)
#' 
#' plot_cluster_network(
#'     object = Seu,
#'     IDC_summary = IDC_summary_scRNA,
#'     edge_df = edge_df
#' )
#' }
#' 
#' # Clustering of scExp scH3K27ac object (Paired-Tag)
#' if(requireNamespace("ChromSCape", quietly=TRUE)){
#' 
#' data("scExp", package = "IDclust")
#' data("IDC_summary_scEpigenomics", package = "IDclust")
#' 
#' plot_cluster_network(
#'     object = scExp,
#'     IDC_summary = IDC_summary_scEpigenomics,
#'     color_by = "IDcluster",
#'     cluster_col = "IDcluster",
#'     colors = NULL,
#'     node_size_factor = 7.5,
#'     edge_size_factor = 1,
#'     function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE,
#'                                                          flip.y = FALSE),
#'     edge_df = topmarkers
#' )
#'
#' # Plotting proportion of cells activating a specific gene in  scExp scH3K27ac 
#' # object (Paired-Tag)
#' plot_cluster_network(
#'     object = scExp,
#'     IDC_summary = IDC_summary_scEpigenomics,
#'     color_by = "Tcf4", # a gene contained in the scExp object
#'     threshold_to_define_feature_active = 1,
#'     gene_col = "Gene",
#'     max_distanceToTSS = 1000,
#'     cluster_col = "IDcluster",
#'     colors = NULL,
#'     node_size_factor = 7.5,
#'     edge_size_factor = 1,
#'     function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE,
#'                                                          flip.y = FALSE)
#' )
#' 
#' # Adding on the edges the 3 top markers of each clusters in scExp H3K27ac
#' # object (Paired-Tag)
#' 
#' # Adding gene information in the IDC_DA
#' data("IDC_DA_scEpigenomics", package = "IDclust")
#' IDC_DA_scEpigenomics = add_gene_to_DA_list(
#'     scExp = scExp, 
#'     IDC_DA = IDC_DA_scEpigenomics
#' )
#' 
#' # Finding the 3 top markers per cluster
#' topmarkers = top_differential_markers(
#'     IDC_DA_scEpigenomics,
#'     top = 3,
#'     gene_col = "Gene",
#'     logFC_col = "logFC",
#'     qvalue_col = "qval",
#'     order_by = "logFC_col",
#'     pseudogene_pattern = "Rik|Vmn|Gm|AW"
#' )
#' 
#' # Concatenate top 3 markers per cluster/cluster_of_origin 
#' topmarkers = topmarkers %>% dplyr::group_by(cluster_of_origin, cluster) %>%
#'  dplyr::summarise(Term = paste(Gene, collapse = " "))
#' 
#' plot_cluster_network(
#'     object = scExp,
#'     IDC_summary = IDC_summary_scEpigenomics,
#'     edge_df = topmarkers
#' )
#' 
#' }
plot_cluster_network <- function(object, ...) {
    UseMethod(generic = 'plot_cluster_network', object = object)
}

#' Plot Iterative Differential Clustering network
#'
#' @param object A SingleCellExperiment object clustered with
#'  [iterative_differential_clustering()].
#' @param IDC_summary Optional. A data.frame of differential
#' analyses summary outputed by [iterative_differential_clustering()] when 
#' saving option is TRUE. Use to determine the width of the edges based on the
#' number of differential features of the given marker.
#' @param color_by A character specifying the column of the SingleCellExperiment
#'  to use for coloring the nodes.
#' @param cluster_col A character specifying the column of the 
#' SingleCellExperiment to use to store the iterative differential clusters. 
#' @param node_size_factor A numeric specifying the size of the nodes.  
#' @param function_layout A function of g for the layout of the graph.
#' @param threshold_to_define_feature_active If color_by is a gene, an integer
#' specifying the threshold above which a gene is considered as active in any 
#' given cell.
#' @param gene_col If color_by is a gene, a character specifying the column in
#' the rowData of the object
#' @param max_distanceToTSS If color_by is a gene, the maximum distance to TSS
#' to consider a gene  linked to a region. Used only if "color_by" is a gene 
#' name.
#' @param legend A logical indicating whether to plot the legend or not.
#' @param edge_df A data.frame containing column 'Term', 'cluster' and 
#' 'cluster_of_origin'. Typically obtained by [top_enriched_pathways()]. 
#' Will add label (Term) on the corresponding edges (cluster_of_origin to 
#' cluster).
#' @param ... Additional parameters passed to the plot function.
#'
#' @return A hierarchical network of cluster assignation:
#' * Size of nodes reflects the number of cells
#' * Width of edges reflects the number of differential features defining a cluster    
#' * Color of nodes reflects the repartition of cells according to 'color_by'
#'
#' @importFrom igraph layout_as_tree simplify graph_from_adjacency_matrix
#' get.edges
#' @importFrom grDevices colors
#' @export
#'
#' @rdname plot_cluster_network
#' @exportS3Method plot_cluster_network default
#' 
plot_cluster_network.default <- function(
    object,
    IDC_summary,
    color_by = "IDcluster",
    cluster_col = "IDcluster",
    colors = NULL,
    node_size_factor = 7.5,
    edge_size_factor = 1,
    threshold_to_define_feature_active = 1,
    max_distanceToTSS = 1000,
    gene_col = "Gene",
    function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE, flip.y = FALSE),
    legend = TRUE,
    edge_df = NULL,
    ...
){
    categ = TRUE
    if(!color_by %in% colnames(SingleCellExperiment::colData(object))) {
        if(color_by %in% unlist(strsplit(SingleCellExperiment::rowData(object)[[gene_col]], split = ", ", fixed = T)))
            categ = FALSE else stop("color_by must be in colnames",
            "(object@meta.data) or in rownames(object).")
    }
    
    object[[cluster_col]] = gsub("Omega:","", object[[cluster_col]])
    if(!is.null(IDC_summary)){
        IDC_summary$true_subcluster = gsub("Omega:","", IDC_summary$true_subcluster)
        IDC_summary$cluster_of_origin = gsub("Omega:","", IDC_summary$cluster_of_origin)
    }

    if(is.null(colors) & categ) {
        colors = c("#4285F4", "#DB4437", "#F4B400", "#0F9D58", "slategray",
                   grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)])
        colors = colors[seq_along(unique(unlist(object[[color_by]])))]
    }
        
    
    annot = as.data.frame(SingleCellExperiment::colData(object))
    
    # Make a color data.frame
    if(categ){
        color_df = data.frame(
            "color_by" = sort(unique(object[[color_by]])),
            "color_by_color" = colors
        )
    } else {
        annot_feature = as.data.frame(SummarizedExperiment::rowRanges(object))
        annot_feature = annot_feature %>% 
            tidyr::separate_rows(.data[[gene_col]], sep = ", ") %>% 
            dplyr::group_by(.data[[gene_col]]) %>%
            dplyr::slice_min(.data[["distanceToTSS"]])
        annot_feature = annot_feature %>%
            dplyr::filter( .data[["distanceToTSS"]] < max_distanceToTSS)
        
        counts = SingleCellExperiment::counts(object[annot_feature$ID[which(
            annot_feature[[gene_col]] == color_by)],])
        
        if(!is.null(counts) & nrow(counts) > 1 ){
            counts = Matrix::rowSums(counts)
        }
        counts = as.numeric(counts)
        sel = which(counts >= threshold_to_define_feature_active)
        counts[] = "Inactive"
        counts[sel] = "Active"
        
        color_df = data.frame(
            "color_by" = sort(unique(counts)),
            "color_by_color" = c("red","grey85")
        )
        annot[, color_by] = counts
    }
    colnames(color_df) = gsub("color_by", color_by, colnames(color_df))
    
    clusters = object[[cluster_col]]
    
    annot$partition_0 = "Omega"
    
    cluster_list = sapply(clusters, function(i) strsplit(i, split = ":", fixed =TRUE))
    max_partition_depth = max(sapply(cluster_list, length))
    for(i in seq_len(max_partition_depth)){
        annot[,paste0("partition_",i)] = sapply(cluster_list, function(l){
            paste(l[1:min(i, length(l)) ], collapse = ":") })
        na_idxs = which(is.na(annot[,paste0("partition_",i)]))
        if(i > 1 & length(na_idxs)) annot[na_idxs, paste0("partition_",i)] =
            annot[na_idxs, paste0("partition_",i-1)]
    }
    
    # Table of cell and cluster assignation for each 'partition depth'
    repartition = annot
    repartition = repartition[,c(paste0("partition_",0:max_partition_depth), color_by)]
    
    # Define all final clusters
    all_nodes = sapply(paste0("partition_",0:max_partition_depth), 
                       function(j){unique(repartition[,j])})
    all_nodes = unique(unlist(all_nodes))
    
    # Create adjency matrix for graph
    adj.mat = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                     dimnames = list(all_nodes, all_nodes))
    adj.mat[1, unique(repartition$partition_1)] = 1
    
    df = data.frame(name =all_nodes)
    rownames(df) = df$name
    
    df$size = 0
    df$ndiff = 0
    df$isFinal = 2
    prop.base_clust = list()
    empty = table(repartition[[color_by]])
    empty[] = 0
    
    # Fill adjency matrix with link between nodes (clusters).
    # Size of nodes reflects the number of cells
    # Width of edges reflects the number of differential features defining a cluster
    # Color of nodes reflects the repartition of cells according to the 'color_by'
    # column
    for(level in 1:(ncol(repartition)-1)){
        for(i in unique(repartition[,level])){
            df$size[df$name == i] = length(which(repartition[,level] == i))
            if(i %in% repartition[,ncol(repartition)-1]) df$isFinal[df$name == i] = 1
            # if(is.null(prop.base_clust[[i]])){
            prop.base_clust[[i]] = empty
            tab = table(repartition[which(repartition[,level] == i), color_by])
            prop.base_clust[[i]][match(names(tab), names(empty))] = as.numeric(tab)
            prop.base_clust[[i]] = as.numeric(prop.base_clust[[i]])
            # } 
            if(i != "Omega") df$ndiff[df$name == i] = IDC_summary$n_differential[which(IDC_summary$true_subcluster == i)][1]
            if(level != (ncol(repartition)-1)) adj.mat[i, unique(repartition[,level + 1 ][which(repartition[,level] == i)])] = 1
        }
        
    }
    
    # Make sure that all the edge have a minimum width of 1
    df$ndiff[which(df$ndiff == 0)] = 1
    df$size = node_size_factor * sqrt(50 * df$size / nrow(repartition))
    df$size[1] = df$size[1] / 1.5
    df$ndiff = edge_size_factor * log2(df$ndiff+1)
    g = igraph::simplify( igraph::graph_from_adjacency_matrix(adjmatrix = adj.mat))
    
    edges_ids =  igraph::get.edges(g, es = 1:(nrow(df)-1))
    
    if(legend) par(mar = c(9, 0, 4, 8) + 0.1)
    
    df$pathway = ""
    if(!is.null(edge_df)){
      edge_df$Term = gsub(" \\(.*", "", edge_df$Term)
      edge_df = edge_df %>% dplyr::group_by(cluster, cluster_of_origin) %>% 
          dplyr::summarise(Term = head(Term, 1))
      df$pathway = edge_df$Term[match(df$name, gsub("Omega:", "",paste0(edge_df$cluster_of_origin, ":", edge_df$cluster)))]
    }
    
    layout = function_layout(g)
    plot(g,
         layout = layout,
         vertex.shape = c("pie"),             # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
         vertex.pie = prop.base_clust,
         rescale=TRUE,
         vertex.pie.lty = df$isFinal,
         vertex.pie.color=list(color_df[,paste0(color_by, "_color")]),
         vertex.size = df$size,                          # Size of the node (default is 15)
         vertex.size2 = NA,
         edge.width=df$ndiff[edges_ids[,2]],                        # Edge width, defaults to 1
         edge.arrow.size=0,                           # Arrow size, defaults to 1
         edge.arrow.width=0,                          # Arrow width, defaults to 1
         edge.lty=c("solid"),
         vertex.label.color=c("black"),
         vertex.label.family='sans',                   # Font family of the label (e.g.“Times”, “Helvetica”)
         vertex.label.font=c(1),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
         vertex.label.cex=0.125 * log2(df$size),                 # Font size (multiplication factor, device-dependent)
         vertex.label.dist=c(3,rep(1, nrow(df) -1)),                           # Distance between the label and the vertex
         vertex.label.degree=c(4.7,rep(90, nrow(df) -1)), 
         margin = c(-0.2,-1,-0.2,-1),
         main = color_by,
         ...)

    arguments <- list(...)
    if(!"edge.label.cex" %in% names(arguments)) edge.label.cex = 0.5 else
      edge.label.cex = arguments$edge.label.cex
    

    plot(g,
         layout = layout,
         vertex.shape="none", 
         edge.arrow.size=0,                           # Arrow size, defaults to 1
         edge.arrow.width=0,
         edge.lty=0,
         edge.width=0,
         vertex.label.cex = 0,
         edge.label = lapply(df$pathway[edges_ids[,2]], function(s) paste(strwrap(s, width=20), collapse = "\n")),
         edge.label.family='sans',     
         edge.label.cex = edge.label.cex,
         vertex.label="",  
         edge.label.color='black', 
         add = TRUE,
         ...)
    
    if(legend){
        par(mar = c(5, 4, 4, 2) + 0.1)
        sizes = c(50, 100, 200, 500, 1000)
        sizes. = node_size_factor * sqrt(50 * sizes / nrow(repartition))
        
        legend(x = 0.5,
               y = -1.25,
               ncol = 6,
               xjust = 0.5,
               box.lty=0,
               legend=unique(color_df[,1]),
               pch=20,
               pt.cex = 4,
               y.intersp = 2,
               col=as.vector(unique(color_df[,2])),
               xpd = TRUE)

        legend(x = 1.5,
               y =1.35,
               title = "#Ncells",
               legend=as.character(sizes),
               box.lty=0,
               pch=1,
               cex = 1,
               pt.cex = sizes./10,
               y.intersp = 1.5,
               x.intersp = 2,
               col="black",
               xpd = TRUE)

        
        text(c(1.5,1.5)  + 0.35, c(0.2,-0.05),
             labels = c("Final", "Transient"))
        
        circle. <- function (r, x0, y0,  ...){
            t <- seq(0, 2 * pi, by = 0.01)
            x <- r * cos(t) + x0
            y <- r * sin(t) + y0
            lines(x, y, ...)
        }
        circle.(r = 0.1,
               x0 =  1.5,
               y0 = 0.2,
               lty = 1, 
               col = "black",
               xpd =TRUE)
        circle.(r = 0.1,
                x0 = 1.5,
                y0 =  -0.05,
                lty = 2, 
                col = "black",
                xpd =TRUE)
        
        if(!is.null(IDC_summary)){
            sizes = c(5,10,20,50,100)
            sizes. = edge_size_factor * log2(sizes+1)
            legend(x = 1.5,
                   y = -0.35,
                   title = "#Nmarkers",
                   legend=as.character(sizes),
                   box.lty=0,
                   lty = 1,
                   cex = 1,
                   lwd = 1.5*sizes.,
                   col = "grey65",
                   y.intersp = c(1.35),
                   x.intersp = 2,
                   xpd = TRUE
                   )
        }
    }
    
    return(color_df)
}


#' Plot Iterative Differential Clustering network
#'
#' @param object A Seurat object clustered with [iterative_differential_clustering()]
#' @param IDC_summary Optional. A data.frame of differential
#' analyses summary outputed by [iterative_differential_clustering()] when 
#' saving option is TRUE. Use to determine the width of the edges based on the
#' number of differential features of the given marker.
#' @param color_by A character specifying the column of the Seurat
#'  to use for coloring the nodes.
#' @param cluster_col A character specifying the column of the 
#' Seurat to use to store the iterative differential clusters. 
#' @param colors A character vector of colors. If NULL, will take R default 
#' color.
#' @param node_size_factor A numeric specifying  a multiplicator of the size of
#' the nodes.  
#' @param edge_size_factor A numeric specifying a multiplicator of the 
#' size of the edges.  
#' @param function_layout A function of g for the layout of the graph.
#' @param threshold_to_define_feature_active If color_by is a gene, an integer
#' specifying the threshold above which a gene is considered as active in any 
#' given cell.
#' @param assay If color_by is a gene, the assay in which to retrieve the counts.
#' @param legend A logical indicating whether to plot the legend or not.
#' @param edge_df (Optional). A data.frame obtained by
#'  [top_enriched_pathways()] containing the top 1 pathway enriched per cluster
#'  to display it on the edges.
#' @param ... Additional parameters passed to the plot function.
#' @return A hierarchical network of cluster assignation:
#' * Size of nodes reflects the number of cells
#' * Width of edges reflects the number of differential features defining a cluster    
#' * Color of nodes reflects the repartition of cells according to 'color_by'
#'
#' @importFrom igraph layout_as_tree simplify graph_from_adjacency_matrix
#' get.edges
#' @importFrom grDevices colors
#' @export
#'
#' @rdname plot_cluster_network
#' @exportS3Method plot_cluster_network Seurat
#' 
plot_cluster_network.Seurat <- function(
    object,
    IDC_summary = NULL,
    color_by = "IDcluster",
    cluster_col = "IDcluster",
    colors = NULL,
    node_size_factor = 7.5,
    edge_size_factor = 1,
    threshold_to_define_feature_active = 2,
    function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE, flip.y = FALSE),
    assay = "RNA",
    legend = TRUE,
    edge_df = NULL,
    ...
){
    categ = TRUE
    if(!color_by %in% colnames(object@meta.data)) {
        if(color_by %in% rownames(object)) categ = FALSE else stop("color_by must ",
                                                                   "be in colnames(object@meta.data) or in rownames(object).")
    }
    
    object@meta.data[,cluster_col] = gsub("Omega:","", object@meta.data[,cluster_col])
    
    if(!is.null(IDC_summary)){
        IDC_summary$true_subcluster = gsub("Omega:","", IDC_summary$true_subcluster)
        IDC_summary$cluster_of_origin = gsub("Omega:","", IDC_summary$cluster_of_origin)
    }

    if(is.null(colors)){
        colors = c("#4285F4", "#DB4437", "#F4B400", "#0F9D58", "slategray",
                   grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)])
        colors = colors[seq_along(unique(unlist(object[[color_by]])))]
    }
    
    # Annotation
    annot = as.data.frame(object@meta.data)
    
    # Make a color data.frame
    if(categ){
        color_df = data.frame(
            "color_by" = sort(unique(object@meta.data[,color_by])),
            "color_by_color" = colors
        )
    } else {
        counts = as.numeric(Seu@assays[[assay]]@counts[color_by,])
        sel = which(counts >= threshold_to_define_feature_active)
        counts[] = "Inactive"
        counts[sel] = "Active"
        
        color_df = data.frame(
            "color_by" = sort(unique(counts)),
            "color_by_color" = c("red","grey85")
        )
        annot[, color_by] = counts
    }
    colnames(color_df) = gsub("color_by", color_by, colnames(color_df))
    
    
    clusters = object[[cluster_col]]
    
    annot$partition_0 = "Omega"
    
    cluster_list = sapply(clusters, function(i) strsplit(i, split = ":", fixed =TRUE))
    max_partition_depth = max(sapply(cluster_list, length))
    for(i in seq_len(max_partition_depth)){
        annot[,paste0("partition_",i)] = sapply(cluster_list, function(l){
            paste(l[1:min(i, length(l)) ], collapse = ":") })
        na_idxs = which(is.na(annot[,paste0("partition_",i)]))
        if(i > 1 & length(na_idxs)) annot[na_idxs, paste0("partition_",i)] =
            annot[na_idxs, paste0("partition_",i-1)]
    }
    
    # Table of cell and cluster assignation for each 'partition depth'
    repartition = annot
    repartition = repartition[,c(paste0("partition_",0:max_partition_depth), color_by)]
    
    # Define all final clusters
    all_nodes = sapply(paste0("partition_",0:max_partition_depth), 
                       function(j){unique(repartition[,j])})
    all_nodes = unique(unlist(all_nodes))
    
    # Create adjency matrix for graph
    adj.mat = matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                     dimnames = list(all_nodes, all_nodes))
    adj.mat[1, unique(repartition$partition_1)] = 1
    
    df = data.frame(name =all_nodes)
    rownames(df) = df$name
    
    df$size = 0
    df$ndiff = 0
    df$isFinal = 2
    prop.base_clust = list()
    empty = table(repartition[[color_by]])
    empty[] = 0
    
    # Fill adjency matrix with link between nodes (clusters).
    # Size of nodes reflects the number of cells
    # Width of edges reflects the number of differential features defining a cluster
    # Color of nodes reflects the repartition of cells according to the 'color_by'
    # column
    for(level in 1:(ncol(repartition)-1)){
        for(i in unique(repartition[,level])){
            df$size[df$name == i] = length(which(repartition[,level] == i))
            if(i %in% repartition[,ncol(repartition)-1]) df$isFinal[df$name == i] = 1
            # if(is.null(prop.base_clust[[i]])){
            prop.base_clust[[i]] = empty
            tab = table(repartition[which(repartition[,level] == i), color_by])
            prop.base_clust[[i]][match(names(tab), names(empty))] = as.numeric(tab)
            prop.base_clust[[i]] = as.numeric(prop.base_clust[[i]])
            # } 
            if(!is.null(IDC_summary)) {
                if(i != "Omega") df$ndiff[df$name == i] = IDC_summary$n_differential[which(IDC_summary$true_subcluster == i)][1]
            } else {
                if(i != "Omega") df$ndiff[df$name == i] = 20
            }
            if(level != (ncol(repartition)-1)) adj.mat[i, unique(repartition[,level + 1 ][which(repartition[,level] == i)])] = 1
        }
        
    }
    # Make sure that all the edge have a minimum width of 1
    df$ndiff[which(df$ndiff == 0)] = 1
    df$size = node_size_factor * sqrt(50 * df$size / nrow(repartition))
    df$size[1] = df$size[1] / 1.5
    df$ndiff = edge_size_factor * log2(df$ndiff+1)
    g = igraph::simplify( igraph::graph_from_adjacency_matrix(adjmatrix = adj.mat))
    
    edges_ids =  igraph::get.edges(g, es = 1:(nrow(df)-1))
    
    if(legend) par(mar = c(9, 0, 4, 8) + 0.1)
    
    df$pathway = ""
    if(!is.null(edge_df)){
        edge_df$Term = gsub(" \\(.*", "", edge_df$Term)
        edge_df = edge_df %>% group_by(cluster, cluster_of_origin) %>% 
            dplyr::summarise(Term = head(Term, 1))
        df$pathway = edge_df$Term[match(df$name, gsub("Omega:", "",paste0(edge_df$cluster_of_origin, ":", edge_df$cluster)))]
    }
    

    layout = function_layout(g)
    plot(g,
         layout = layout,
         vertex.shape = c("pie"),             # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
         vertex.pie = prop.base_clust,
         rescale=TRUE,
         vertex.pie.lty = df$isFinal,
         vertex.pie.color=list(color_df[,paste0(color_by, "_color")]),
         vertex.size = df$size,                          # Size of the node (default is 15)
         vertex.size2 = NA,
         edge.width=df$ndiff[edges_ids[,2]],                        # Edge width, defaults to 1
         edge.arrow.size=0,                           # Arrow size, defaults to 1
         edge.arrow.width=0,                          # Arrow width, defaults to 1
         edge.lty=c("solid"),
         vertex.label.color=c("black"),
         vertex.label.family='sans',                   # Font family of the label (e.g.“Times”, “Helvetica”)
         vertex.label.font=c(1),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
         vertex.label.cex=0.125 * log2(df$size),                 # Font size (multiplication factor, device-dependent)
         vertex.label.dist=c(3,rep(1, nrow(df) -1)),                           # Distance between the label and the vertex
         vertex.label.degree=c(4.7,rep(90, nrow(df) -1)), 
         margin = c(-0.2,-1,-0.2,-1),
         main = color_by,
         ...)
    
    arguments <- list(...)
    if(!"edge.label.cex" %in% names(arguments)) edge.label.cex = 0.5 else
      edge.label.cex = arguments$edge.label.cex
    
    plot(g,
         layout = layout,
         vertex.shape="none", 
         edge.arrow.size=0,                           # Arrow size, defaults to 1
         edge.arrow.width=0,
         edge.lty=0,
         edge.width=0,
         vertex.label.cex = 0,
         edge.label = lapply(df$pathway[edges_ids[,2]], function(s) paste(strwrap(s, width=20), collapse = "\n")),
         edge.label.family='sans',     
         vertex.label="",  
         edge.label.color='black', 
         add = TRUE,
         edge.label.cex = edge.label.cex,
         ...)

    if(legend){
        par(mar = c(5, 4, 4, 2) + 0.1)
        sizes = c(50, 100, 200, 500, 1000)
        sizes. = node_size_factor * sqrt(50 * sizes / nrow(repartition))
        
        legend(x = 0.5,
               y = -1.25,
               ncol = 6,
               box.lty=0,
               legend=unique(color_df[,1]),
               pch=20,
               xjust = 0.5,
               pt.cex = 4,
               y.intersp = 2,
               col=as.vector(unique(color_df[,2])),
               xpd = TRUE)
        
        legend(x = 1.5,
               y =1.35,
               title = "#Ncells",
               legend=as.character(sizes),
               box.lty=0,
               pch=1,
               cex = 1,
               pt.cex = sizes./10,
               y.intersp = 1.5,
               x.intersp = 2,
               col="black",
               xpd = TRUE)
        
        
        text(c(1.5,1.5)  + 0.35, c(0.2,-0.05),
             labels = c("Final", "Transient"))
        
        circle. <- function (r, x0, y0,  ...){
            t <- seq(0, 2 * pi, by = 0.01)
            x <- r * cos(t) + x0
            y <- r * sin(t) + y0
            lines(x, y, ...)
        }
        circle.(r = 0.1,
                x0 =  1.5,
                y0 = 0.2,
                lty = 1, 
                col = "black",
                xpd =TRUE)
        circle.(r = 0.1,
                x0 = 1.5,
                y0 =  -0.05,
                lty = 2, 
                col = "black",
                xpd =TRUE)
        
        if(!is.null(IDC_summary)){
            sizes = c(5,10,20,50,100)
            sizes. = edge_size_factor * log2(sizes+1)
            legend(x = 1.5,
                   y = -0.35,
                   title = "#Nmarkers",
                   legend=as.character(sizes),
                   box.lty=0,
                   lty = 1,
                   cex = 1,
                   lwd = 1.5*sizes.,
                   col = "grey65",
                   y.intersp = c(1.35),
                   x.intersp = 2,
                   xpd = TRUE
            )
        }
    }
    
    return(color_df)
}

