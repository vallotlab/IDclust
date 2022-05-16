#' Plot Iterative Differential Clustering network
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_cluster_network <- function(object, ...) {
    UseMethod(generic = 'plot_cluster_network', object = object)
}

#' Plot Iterative Differential Clustering network
#'
#' @param object A SingleCellExperiment object clustered with [iterative_differential_clustering()]
#' @param differential_summary_df Optional. A data.frame of differential
#' analyses summary outputed by [iterative_differential_clustering()] when 
#' saving option is TRUE. Use to determine the width of the edges based on the
#' number of differential features of the given marker.
#' @param color_by A character specifying the column of the SingleCellExperiment
#'  to use for coloring the nodes.
#' @param cluster_col A character specifying the column of the 
#' SingleCellExperiment to use to store the iterative differential clusters. 
#' @param node_size_factor A numeric specifying the size of the nodes.  
#' @param function_layout A function of g for the layout of the graph.
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
#' @method plot_cluster_network default
#' @S3method plot_cluster_network default
#' @examples
plot_cluster_network.default <- function(
    object,
    differential_summary_df,
    color_by = "cell_cluster",
    cluster_col = "cell_cluster",
    colors = NULL,
    node_size_factor = 7.5,
    edge_size_factor = 1,
    function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE, flip.y = F),
    ...
){
    
  object[[cluster_col]] = gsub("Omega:","", object[[cluster_col]])
  differential_summary_df$true_subcluster = gsub("Omega:","", differential_summary_df$true_subcluster)
  differential_summary_df$cluster_of_origin = gsub("Omega:","", differential_summary_df$cluster_of_origin)
  if(is.null(colors)) 
    colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
  
  # Make a color data.frame
  color_df = data.frame(
    "color_by" = sort(unique(object[[color_by]])),
    "color_by_color" = sample(colors, length(unique(object[[color_by]])), replace = F)
  )
  colnames(color_df) = gsub("color_by", color_by, colnames(color_df))
  
  annot = as.data.frame(SingleCellExperiment::colData(object))
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
  empty = table(object[[color_by]])
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
      if(i != "Omega") df$ndiff[df$name == i] = differential_summary_df$n_differential[which(differential_summary_df$true_subcluster == i)][1]
      if(level != (ncol(repartition)-1)) adj.mat[i, unique(repartition[,level + 1 ][which(repartition[,level] == i)])] = 1
    }
    
  }
  # Make sure that all the edge have a minimum width of 1
  df$ndiff[which(df$ndiff == 0)] = 1
  df$size = node_size_factor * sqrt(50 * df$size / nrow(repartition))
  df$ndiff = edge_size_factor * log2(df$ndiff+1)
  g = igraph::simplify( igraph::graph_from_adjacency_matrix(adjmatrix = adj.mat))
  
  edges_ids =  igraph::get.edges(g, es = 1:(nrow(df)-1))

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
       vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
       vertex.label.font=c(1),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
       vertex.label.cex=0.125 * log2(df$size),                 # Font size (multiplication factor, device-dependent)
       vertex.label.dist=c(3,rep(1, nrow(df) -1)),                           # Distance between the label and the vertex
       vertex.label.degree=c(4.7,rep(90, nrow(df) -1)), 
       margin = c(-0.2,-1,-0.2,-1),
       ...
  )


return(color_df)
}


#' Plot Iterative Differential Clustering network
#'
#' @param object A Seurat object clustered with [iterative_differential_clustering()]
#' @param differential_summary_df Optional. A data.frame of differential
#' analyses summary outputed by [iterative_differential_clustering()] when 
#' saving option is TRUE. Use to determine the width of the edges based on the
#' number of differential features of the given marker.
#' @param color_by A character specifying the column of the Seurat
#'  to use for coloring the nodes.
#' @param cluster_col A character specifying the column of the 
#' Seurat to use to store the iterative differential clusters. 
#' @param node_size_factor A numeric specifying the size of the nodes.  
#' @param function_layout A function of g for the layout of the graph.
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
#' @method plot_cluster_network Seurat
#' @S3method plot_cluster_network Seurat
#' @examples
plot_cluster_network.Seurat <- function(
    object,
    differential_summary_df,
    color_by = "cell_cluster",
    cluster_col = "cell_cluster",
    colors = NULL,
    node_size_factor = 7.5,
    edge_size_factor = 1,
    function_layout = function(g) igraph::layout_as_tree(g, root = 1, circular = TRUE, flip.y = F),
    ...
){
    
    object@meta.data[,cluster_col] = gsub("Omega:","", object@meta.data[,cluster_col])
    differential_summary_df$true_subcluster = gsub("Omega:","", differential_summary_df$true_subcluster)
    differential_summary_df$cluster_of_origin = gsub("Omega:","", differential_summary_df$cluster_of_origin)
    if(is.null(colors)) 
        colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
    
    # Make a color data.frame
    color_df = data.frame(
        "color_by" = sort(unique(object@meta.data[,color_by])),
        "color_by_color" = sample(colors, length(unique(object@meta.data[,color_by])), replace = F)
    )
    colnames(color_df) = gsub("color_by", color_by, colnames(color_df))
    
    annot = as.data.frame(object@meta.data)
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
    empty = table(object[[color_by]])
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
            if(i != "Omega") df$ndiff[df$name == i] = differential_summary_df$n_differential[which(differential_summary_df$true_subcluster == i)][1]
            if(level != (ncol(repartition)-1)) adj.mat[i, unique(repartition[,level + 1 ][which(repartition[,level] == i)])] = 1
        }
        
    }
    # Make sure that all the edge have a minimum width of 1
    df$ndiff[which(df$ndiff == 0)] = 1
    df$size = node_size_factor * sqrt(50 * df$size / nrow(repartition))
    df$ndiff = edge_size_factor * log2(df$ndiff+1)
    g = igraph::simplify( igraph::graph_from_adjacency_matrix(adjmatrix = adj.mat))
    
    edges_ids =  igraph::get.edges(g, es = 1:(nrow(df)-1))
    
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
         vertex.label.family="Helvetica",                   # Font family of the label (e.g.“Times”, “Helvetica”)
         vertex.label.font=c(1),                  # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
         vertex.label.cex=0.125 * log2(df$size),                 # Font size (multiplication factor, device-dependent)
         vertex.label.dist=c(3,rep(1, nrow(df) -1)),                           # Distance between the label and the vertex
         vertex.label.degree=c(4.7,rep(90, nrow(df) -1)), 
         margin = c(-0.2,-1,-0.2,-1),
         ...
    )
    
    
    return(color_df)
}

