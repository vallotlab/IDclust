args <- commandArgs(TRUE)
input = list()
input$scExp_path = as.character(args[1])
input$working_directory = as.character(args[2])
input$output_dir = as.character(args[3])
input$parameter_number = as.numeric(args[4])

# Directories -------------------------------------------------------------
maindir= input$working_directory
source(file.path(maindir, "scripts", "global_var_cluster.R"))
source(file.path(maindir, "scripts", "functions.R"))
source(file.path(maindir, "scripts", "differential_clustering_scHistone.R"))

scExp = qs::qread(input$scExp_path)


parameter_set = read.csv(file.path(maindir, "output", "Marks", "PairedTag", "differential_clustering_scHistone_parameters_optimal_parameters.csv"))
limit_by_proportion = qs::qread(file.path(maindir, "output", "Marks", "PairedTag", parameter_set$mark[input$parameter_number], "FDR_proportion.qs"))

cat("Doing ",  basename(input$scExp_path), "with\n",
    parameter_set$name[input$parameter_number], "parameters:\n",
    "- bins =", parameter_set$bins[input$parameter_number], "(nrow =",nrow(scExp),")\n",
    "- starting_res =", parameter_set$start.res[input$parameter_number], "\n",
    "- nPCA =", parameter_set$nPCA[input$parameter_number],  "\n",
    "- percent_feature =", parameter_set$percent.feature[input$parameter_number], "\n",
    "- qval.th =", parameter_set$qval.th[input$parameter_number], "\n",
    "- min.pct.cell_assigned =", parameter_set$min.pct.cell_assigned[input$parameter_number], "\n",
    "- quantile.activation =", parameter_set$quantile.activation[input$parameter_number], "\n",
    "- FC =", parameter_set$FC[input$parameter_number], "\n"
    )


scExp$cell_cluster = gsub("C", "A", as.character(find_clusters_louvain_Seurat_scExp(scExp,
                                                                                    k = 100,
                                                                                    nPCA = 10 - 1,
                                                                                    resolution = parameter_set$start.res[input$parameter_number],
                                                                                    use.dimred = "PCA")))
starting_differential_summary_df = data.frame("cluster" = unique(scExp$cell_cluster), "n_differential" = 10) # we assume the initial clusters are differential

# Run differential clustering
scExp_differential_clustering(scExp = scExp,
                              starting_differential_summary_df = starting_differential_summary_df,
                              output_dir = input$output_dir,
                              nPCA = parameter_set$nPCA[input$parameter_number],
                              percent_feature = parameter_set$percent.feature[input$parameter_number],
                              quantile.activation = parameter_set$quantile.activation[input$parameter_number],
                              FC.th = parameter_set$FC[input$parameter_number],
                              qval.th = parameter_set$qval.th[input$parameter_number],
                              min.pct.cell_assigned =  parameter_set$min.pct.cell_assigned[input$parameter_number],
                              limit = parameter_set$limit[input$parameter_number],
                              limit_by_proportion = limit_by_proportion,
                              k = parameter_set$k[input$parameter_number],
                              resolution = parameter_set$resolution[input$parameter_number],
                              runFDR = as.logical(parameter_set$runFDR[input$parameter_number]))
gc()