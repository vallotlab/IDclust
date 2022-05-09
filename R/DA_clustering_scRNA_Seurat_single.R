args <- commandArgs(TRUE)
input = list()
input$Seu_path = as.character(args[1])
input$working_directory = as.character(args[2])
input$output_dir = as.character(args[3])
input$parameter_number = as.numeric(args[4])

# Directories -------------------------------------------------------------
maindir= input$working_directory
source(file.path(maindir, "scripts", "global_var_cluster.R"))
source(file.path(maindir, "scripts", "functions.R"))
source(file.path(maindir, "scripts", "differential_clustering_scRNA_Seurat.R"))

Seu = qs::qread(input$Seu_path)

parameter_set = read.csv(file.path(maindir, "output", "Marks", "PairedTag", "differential_clustering_scRNA_parameters_Seurat.csv"))

cat("Doing ",  basename(input$Seu_path), "with\n",
    parameter_set$name[input$parameter_number], "parameters:\n",
    "- starting_res =", parameter_set$start.res[input$parameter_number], "\n",
    "- nPCA =", parameter_set$nPCA[input$parameter_number],  "\n",
    "- nfeatures =", parameter_set$nfeatures[input$parameter_number], "\n",
    "- logFC.th =",parameter_set$logFC.th[input$parameter_number], "\n",
    "- qval.th =", parameter_set$qval.th[input$parameter_number], "\n",
    "- min.pct =", parameter_set$min.pct[input$parameter_number], "\n",
    "- min.pct.cell_assigned =", parameter_set$min.pct.cell_assigned[input$parameter_number], "\n")

Seu$seurat_clusters = find_clusters_louvain_Seurat(Seu,
                                                   k = parameter_set$start.k[input$parameter_number],
                                                   resolution = parameter_set$start.res[input$parameter_number],
                                                   use.dimred = "pca")

Seu$seurat_clusters <- paste0("A",as.numeric(Seu$seurat_clusters))
Idents(Seu) = Seu$seurat_clusters

starting_differential_summary_df = data.frame("cluster" = unique(Seu$seurat_clusters), "n_differential" = parameter_set$limit[input$parameter_number]) # we assume the initial clusters are differential

print(starting_differential_summary_df)

# Run differential clustering
Seurat_differential_clustering(Seu = Seu,
                               starting_differential_summary_df = starting_differential_summary_df,
                               output_dir = input$output_dir,
                               nPCA = parameter_set$nPCA[input$parameter_number],
                               nfeatures = parameter_set$nfeatures[input$parameter_number],
                               logFC.th = parameter_set$logFC.th[input$parameter_number],
                               qval.th = parameter_set$qval.th[input$parameter_number],
                               min.pct = parameter_set$min.pct[input$parameter_number],
                               min.pct.cell_assigned = parameter_set$min.pct.cell_assigned[input$parameter_number],
                               limit = parameter_set$limit[input$parameter_number],
                               k = parameter_set$k[input$parameter_number],
                               resolution = parameter_set$resolution[input$parameter_number],
                               test.use = "edgeR-LRT",
                               biological_replicate_col = "Target")
