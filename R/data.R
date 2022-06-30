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

#' A data.frame containing the differential markers for scRNA dataset
#' (Paired-Tag - scRNA) 
#' @description A data.frame of differential markers for scRNA
#'
#' @usage data("IDC_DA_scRNA")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format IDC_DA_scRNA - A data.frame obtained by 
#' [iterative_differential_clustering()].
#' 
"IDC_DA_scRNA"


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


#' A data.frame containing the differential markers for scEpigenomics dataset
#' (Paired-Tag - scEpigenomics) 
#' @description A data.frame of differential markers for scEpigenomics
#'
#' @usage data("IDC_DA_scEpigenomics")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format IDC_DA_scEpigenomics - A data.frame obtained by 
#' [iterative_differential_clustering()].
#' 
"IDC_DA_scEpigenomics"


#' A list of linear models of false positive differential markers accross number
#' of cells
#'
#' 
#' @description A list of linear mode#' according to number of cells for scEpigenomics dataset accross various 
#' binsizes (200000, 100000, 50000, 20000, 10000, 5000, TSS, genebody_TSS)
#' and was calculated aggregating false positive values accross multiple
#' marks (H3K9me3, H3K27me3, H3K27ac, H3K4me1, H3K4me3.)
#'
#' @usage data("lm_list")
#' 
#' @references  
#' Joint profiling of histone modifications and transcriptome in single cells
#' from mouse brain, Chenxu Zhu, Yanxiao Zhang,  Yang Eric Li, Jacinta Lucero,
#' . Margarita Behrens, Bing Ren, Nature Methods, 2021
#' \href{https://doi.org/10.1038/s41592-021-01060-3}{Paired-Tag}
#' 
#' @format lm_list - A list of linear models. See [stats::lm()].
#' 
"lm_list"