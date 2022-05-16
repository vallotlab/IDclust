context("Testing IDclust scEpigenomics")

# Functions for testing purposes
if(requireNamespace("ChromSCape")){
    set.seed(47)
    out = create_scDataset_raw(featureType = "window",sparse = TRUE, 
                               batch_id = factor(c(1,1,2,2)))
    mat = out$mat
    annot = out$annot
    batches = out$batches
    
    # Download, extract & format PairedTag dataset - H3K27ac (Zhu et al., 2021)
    temp = tempfile()
    tempdir_1 = tempdir()
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152020&format=file&file=GSE152020%5FPaired%2DTag%5FH3K27ac%5FDNA%5Ffiltered%5Fmatrix%2Etar%2Egz", temp)
    untar(temp, exdir = tempdir_1)
    features = read.table(file.path(tempdir_1, "04.Paired-Tag_H3K27ac_DNA_filtered_matrix", "bins.tsv"),
                          row.names = NULL, header = F, sep = "\t")[,1, drop = F]
    write.table(features, file = file.path(tempdir_1, "04.Paired-Tag_H3K27ac_DNA_filtered_matrix", "features.tsv"),
                row.names = F, col.names = F, quote = F)
    # Reading the matrix with ChromSCape
    out = ChromSCape::read_sparse_matrix(file.path(tempdir_1, "04.Paired-Tag_H3K27ac_DNA_filtered_matrix"),
                                         ref = "mm10", verbose = TRUE)
    unlink(file.path(tempdir_1, "04.Paired-Tag_H3K27ac_DNA_filtered_matrix"), recursive = TRUE)
    
    scExp = ChromSCape::preprocessing_filtering_and_reduction(
        datamatrix = out$datamatrix,
        annot_raw = out$annot_raw,
        min_reads_per_cell = 200,
        max_quantile_read_per_cell = 99,
        n_top_features = nrow(out$datamatrix),
        norm_type = "TFIDF",
        remove_PC = "Component_1",
        subsample_n = NULL,
        ref_genome = "mm10",
        exclude_regions = NULL,
        doBatchCorr = FALSE,
        batch_sels = NULL
        )
    
    scExp = find_clusters_louvain_scExp(scExp, k = 100, resolution = 0.1, use.dimred = "PCA")
    
    outdir = tempdir()
    
    scExp_IDC = iterative_differential_clustering_scEpigenomics(scExp, output_dir = outdir, nPCA = 10, runFDR = F)
    
    #test sparse matrix
    test_that("Sparse matrices", {
        scExp = create_scExp(mat,annot)
        expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
        scExp = filter_scExp(scExp)
        expect_is(SingleCellExperiment::counts(scExp),"dgCMatrix")
        scExp = normalize_scExp(scExp,type = "CPM")
        expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
        scExp=feature_annotation_scExp(scExp)
        expect_is(SummarizedExperiment::rowRanges(scExp),"GRanges")
        scExp = reduce_dims_scExp(scExp,n = 50,batch_correction = FALSE)
        expect_is(SingleCellExperiment::reducedDim(scExp,"PCA"),"data.frame")
        
        scExp = colors_scExp(scExp,annotCol = c("sample_id","batch_id","total_counts"))
        plot_reduced_dim_scExp(scExp,reduced_dim = "PCA",color_by = "sample_id")
        
        scExp = correlation_and_hierarchical_clust_scExp(scExp)
        expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
        scExp = filter_correlated_cell_scExp(scExp)
        expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
        scExp = consensus_clustering_scExp(scExp)
        expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
        expect_is(scExp@metadata$consclust,"list")
        expect_is(scExp@metadata$consclust[[2]]$consensusClass,"integer")
        scExp = choose_cluster_scExp(scExp,nclust = 2)
        expect_is(SingleCellExperiment::normcounts(scExp),"dgCMatrix")
    })
} else{
    message("Testing IDclust scEpigenomic - no package ChromSCape. Skipping tests")
}




