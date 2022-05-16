#' A Seurat object containing the scRNA part of Paired-Tag dataset. (Paired-Tag - scRNA) 
#' 
#' @description A Seurat object of a  subset of 2000 x 2000 cells by genes of
#' the \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag} multiomic
#' paper. The metadata contains annotation from the author as well as 
#' iterative differential clusters run on the subset of cells and genes. The
#' counts were analyzed with Seurat package.
#'
#' @usage data("Seu")
#' 
#' @references  
#' 
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @import SeuratObject
#' @format Seu - A Seurat object. See [Seurat::Seurat()].
#' 
"Seu"

#' A data.frame containing the summary of differential markers for scRNA dataset
#' (Paired-Tag - scRNA) 
#' @description A summary of differential markers for scRNA
#'
#' @usage data("IDC_summary_scRNA")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format IDC_summary_scRNA - A data.frame obtained by 
#' [iterative_differential_clustering()].
#' 
"IDC_summary_scRNA"


#' A SingleCellExperiment object containing the scH3K27ac part of Paired-Tag 
#' dataset. (Paired-Tag - H3K27ac) 
#' 
#' @description A SingleCellExperiment object of a  subset of 2000 x 5000 cells 
#' by features of the
#'  \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag} multiomic
#' paper. The metadata contains annotation from the author as well as 
#' iterative differential clusters run on the subset of cells and features The
#' counts were analyzed with ChromSCape package.
#'
#' @usage data("scExp")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format scExp - A SingleCellExperiment object. See 
#' [SingleCellExperiment::SingleCellExperiment].
#' 
"scExp"

#' A data.frame containing the summary of differential markers for scEpigenomics dataset
#' (Paired-Tag - H3K27ac) 
#' @description A summary of differential markers for scEpigenomics (Paired-Tag - H3K27ac) 
#'
#' @usage data("IDC_summary_scEpigenomics")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format IDC_summary_scEpigenomics - A data.frame obtained by 
#' [iterative_differential_clustering()].
#' 
"IDC_summary_scEpigenomics"