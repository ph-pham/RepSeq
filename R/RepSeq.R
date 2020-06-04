#' RepSeq
#'
#' Tools for analyse of HTS RepSeq data.
#'
#' This package is for analysing high-thoughput sequencing repertoire data. 
#' It focus on clonotype data sets pre-processd by ClonotypeR, rTCR and MiXCR. 
#' Clonotype tables from different samples are merged, filtered, normalized and 
#' tested for differential expression between groups of samples. 
#' Diversity indices are also studied.
#' Results could be visualized in different ways.
#
# Imports
#
# @importClassesFrom data.table data.table
#' @import data.table
#' @import utils
#' @import graphics
#' @import pbapply
#' @import parallel
#' @import methods
# @importFrom DESeq2 DESeqDataSetFromMatrix
# @importFrom pheatmap pheatmap
# @importFrom stats as.formula
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom stats median sd as.formula
#' @importFrom grDevices colorRampPalette
# @importFrom stats sd
# @importFrom ade4 s.class
# @importFrom naturalsort naturalsort
# @importFrom scales percent
# @importFrom limma vennDiagram
# @importFrom limma vennCounts
# @importFrom vegan vegdist
#' @useDynLib RepSeq
#' @name RepSeq
NULL


#"_PACKAGE"
#> [1] "_PACKAGE"