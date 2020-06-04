#------------------- plot functions --------------------#
#' plot a heatmap of V-J count 
#'
#' function plots a heatmap of V-J counts for sampleName 
#'
#' @param x an object of class RepSeqExperiment 
#' @param sampleName sample to plot, sample name must be existed among x column names.
#' @param scale type of data to plot: raw count, percentages or count per million.
#' @return a heatmap
#' @export
#' @seealso \code{\link[pheatmap]{pheatmap}}
# @example
# plotCountVJ()
plotCountVJ <- function(x, sampleName = NULL, scale = c("counts", "percent", "cpm")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.") 
    sNames <- rownames(sData(x))
    if (is.null(sampleName) || sampleName == "") {
        index <- sNames[1]
        cat("Plot for the first sample in x:", index, ".\n")
    } else {
        index <- sampleName 
        if (length(sampleName) > 1) {
            index <- sampleName[1]
            cat("Only the first sample is used, sample:", index, ".\n")           
        }
        if (is.na(match(index, sNames))) stop("Sample ", index, " not found in x.")                 
    }
    cts <- data.table::copy(assay(x))
    tmp <- data.table::dcast(cts[lib == index, ], J~V, value.var = "count", fun.aggregate = sum)
    data2plot <- data.frame(tmp, row.names = 1)    
    if(scale == "percent"){
      data2plot <- prop.table(data2plot)
    }
    if(scale == "cpm"){
      data2plot <- prop.table(data2plot)*10^6
    }   
    # print(data2plot)
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        p <- pheatmap::pheatmap(data2plot, main = paste0("sample: ", index), cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 20, silent = T)
    }
    return(p)
}

# plot usage
#
# plots a heatmap of percentages for each feature: VpJ/V/J/V-J of several repertoires.
#
# @param x a matrix of counts features in rows and samples in columns.
# @param graph.title a string . Default is "Segment usage". 
# @param cluster if TRUE, hierarchical clustering tree is show in the sample dimension using Euclidean distance and Ward's criterion for aggregation.
# @param ... others options will be passed through pheatmap
# @return a heatmap
# @seealso \code{\link[pheatmap]{pheatmap}}
#plotUsage <- function(x, graph.title="Segment usage", cluster=TRUE,...) {
#    if (missing(x)) stop("Dat set is missing")
#    if (!is.matrix(x)) x <- as.matrix(x)
#    if (cluster) {
#         pheatmap::pheatmap(x, main=graph.title, cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", clustering_distance_cols = "euclidean",...)
#            } else {
#            pheatmap::pheatmap(x, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, heigh = 3,...)
#        }       
#}

# plot a heatmap of diversity indices.
#
# plotDiversity plots a heatmap of V-J counts for sampleName .
#
# @param x a matrix of count data.
# @param graph.title a a title to graph.
# @param cluster if TRUE, hierarchical clustering tree is show in the sample dimension using Euclidean distance and Ward's criterion for aggregation.
# @param sample.annotation a data frame containing sample-related information, sample name must be existed among x column names.
# @return a heatmap of count
# @export
# @seealso \code{\link[pheatmap]{pheatmap}}
# @example
# plotCountVJ()
#plotDiversity <- function(x, graph.title="Shannon diversity", cluster=TRUE, sample.annotation=NA) {
#    if (missing(x)) stop("Diversity index data set is missing.")
#    if (!is.matrix(x)) x <- as.matrix(x)
#    if (!is.na(sample.annotation)) sample.annotation <- data.frame(sample.annotation) 
#    if (cluster) {
#    pheatmap::pheatmap(x, main=graph.title, 
#        cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D", clustering_distance_cols = "euclidean",
#        annotation_col=sample.annotation)
#    } else {
#        pheatmap::pheatmap(x, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col=sample.annotation)
#    } 
#}

#' plot spectratype
#'
#' function plots spectratype for V genes of one repertoire (distribution of lengths of CDR3 in amio-acid).
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param sampleName sample to plot, if NULL the first sample of x is plotted.
#' @param scale the type of bars in term of counts, percentages or count per million.
#' @return a barplot
#' @export
# @example
# plotSpectratype
plotSpectratype <- function(x, sampleName = NULL, scale = c("counts", "percent", "cpm")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
    CDR3length=pep=lib=CDR3aa=N=percent=cpm <- NULL
    sNames <- rownames(sData(x))
    if (is.null(sampleName)) {
        index <- sNames[1]
        cat("Plot for the first sample in x:", index, ".\n")
    } else {
        index <- sampleName 
        if (length(sampleName) > 1) {
            index <- sampleName[1]
            cat("Only the first sample is used, sample:", index, ".\n")           
            }                  
        if (is.na(match(index, sNames))) {
            stop("Sample ", index, " not found in x.")                 
        }
    }
    cts <- data.table::copy(assay(x))
    # subset count
    scl <- match.arg(scale)
    if (scl == "counts") {
      data2plot <- cts[lib == index, ][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)]#geom_bar
      p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N, fill = V))
    }
    if (scl == "percent") {
      data2plot <- cts[lib == index, ][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N)]#geom_bar
      p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = percent, fill = V)) + 
                ggplot2::scale_y_continuous(labels = scales::percent)
                
    }
    if (scl == "cpm") {
      data2plot <- cts[lib == index, ][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,cpm := prop.table(N)*10^6]#geom_bar
      p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = cpm, fill = V))
    }
    #ggplot
    p <- p + ggplot2::geom_bar(stat = "identity", position="stack", colour = "black") +
             ggplot2::labs(title = paste(index,": V Distribution by aa length"), x = "CDR3 length (aa)", y = scl) +
             ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size=14),
                axis.text = ggplot2::element_text(size=14),
                text = ggplot2::element_text(size=14), 
                axis.text.x = ggplot2::element_text(size=14),
                axis.text.y = ggplot2::element_text(size=14))
    return(p)
}

#' plot Renyi's profiles
#'
#' function plots Renyi's entropy profile for all samples in a RepSeqExperiment object.
#'
#' @param x an object of class RepSeqExperiment
#' @param alpha a vector of parameters for estimating Renyi's entropy.  
#' @param level level of repertoire to estimated 
#' @param colorBy name of a factor in sample information data (run sData(x)) 
#' @return a graph
#' @export
# @example
# plotRenyiProfiles()
plotRenyiProfiles <- function(x, alpha = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level = c("VpJ", "V", "J", "VJ", "CDR3aa"), colorBy = NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (length(alpha) < 2) stop("At least 2 values of alpha is needed.")
    variable=value <- NULL
    sdata <- sData(x)
    sNames <- rownames(sdata) 
    levelChoice <- match.arg(level)
    # compute Renyi
    tmp <- renyiProfiles(x, scales = alpha, level=levelChoice)
    # plot 
    if (is.null(colorBy)) {
      colNames <- colnames(tmp)
      #setnames(data2plot, colnames(data2plot), c("alpha", sNames))
      data2plot <- data.table::melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "lib")
      data2plot[, alpha := as.numeric(as.character(variable))]
      p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = alpha, y = as.numeric(value), color = lib, group = lib)) +
                ggplot2::geom_line(size = 1.5) +
                ggplot2::labs(title = paste("Level ",levelChoice, ": Renyi's Entropy"), y = "Renyi's Entropy") +
                ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
                    text = ggplot2::element_text(size=14),
                    axis.text = ggplot2::element_text(size=14),
                    axis.text.x = ggplot2::element_text(size=14),
                    axis.text.y = ggplot2::element_text(size=14))
        } else {
            if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
            #ggplot           
            data2plot <- data.table::melt(data = tmp, id.vars = "variable", measure.vars = sNames, variable.name = "lib")
            data2plot[, paste(colorBy) := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "lib"]
            data2plot[, alpha := as.numeric(as.character(variable))]
            p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x = "alpha", y = "value", colour = paste(colorBy))) +
                    ggplot2::geom_line(ggplot2::aes(group = lib), size = 1.5) +
                    ggplot2::labs(title = paste("Level ",levelChoice, ": Renyi's Entropy"), y = "Renyi's Entropy") +
                    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
                        text = ggplot2::element_text(size=14), 
                        axis.text = ggplot2::element_text(size=14),
                        axis.text.x = ggplot2::element_text(size=14),
                        axis.text.y = ggplot2::element_text(size=14))
        }
    return(p)
}

#' plot V or J proportion
#'
#' function plots V or J proportion in all samples in a RepSeqExperiment object.
#'
#' @param x an object of class RepSeqExperiment.
#' @param level level of repertoire to be used.
#' @param sampleName the name of sample or lib to be plotted, sample name must be existed among x column names. If NULL the first sample name is used.
#' @param scale type of data to be plotted, counts, percentages or counts per million. 
#' @return a barplot
#' @export
# @example
# plotPropVJ
plotPropVJ <- function(x, level = c("V", "J"), sampleName = NULL, scale = c("counts", "percent", "cpm")) {
    percent=cpm <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")  
    if (is.null(sampleName)) sName <- rownames(sData(x))[1] else sName <- sampleName
    # match argument
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    data2plot <- cts[lib == sName, lapply(.SD, sum), by = levelChoice, .SDcols = "count"]
    if (scale == "counts") {   
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x = levelChoice, y = "count", fill = levelChoice)) 
    }
    if (scale == "percent") {
        data2plot[, percent := prop.table(count)]
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x = levelChoice, y = "percent", fill = levelChoice)) +
             ggplot2::scale_y_continuous(labels = scales::percent)
    }
    if (scale == "cpm") {
        data2plot[, cpm := prop.table(count)*10^6]
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes_string(x = levelChoice, y = "cpm", fill = levelChoice)) 
    }
    p <- p + ggplot2::geom_bar(width = 0.7, stat = "identity", show.legend=F) +
            ggplot2::scale_x_discrete(limits = data2plot[order(-count)][[levelChoice]]) +
            ggplot2::labs(title = paste(sName, ": ", levelChoice, " Distribution")) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
                            axis.text = ggplot2::element_text(size=14),
                            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
                            axis.text.y = ggplot2::element_text(size = 14))
  return(p)
}

#' plot frequency
#'
#' function plot the histogram of counts for a sample broken into 6 intervals 1, ]1, 10], ]10, 100], ]100, 1000], ]1000, 10000], ]10000, more]. 
#' @param x an object of class RepSeqExperiment
#' @param sampleName sample to be plotted. If NULL the first sample name is used.
#' @return a bar plot
#' @export
# @example
plotFreqVpJ <- function(x, sampleName = NULL){
    interval=percent <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (is.null(sampleName)) {
        sName <- rownames(sData(x))[1]
        cat("Plot for the first sample in x:", sName, ".\n")
    } else sName <- sampleName
    data2plot <- data.table::copy(assay(x))
    #print(data2plot)
    f <- function(x){
        if(x == 1) "1"
        else if(x <= 10) "]1, 10]"
        else if(x <= 100) "]10, 100]"
        else if(x <= 1000) "]100, 1000]"
        else if(x <= 10000) "]1000, 10000]"
        else "]10000, 100000]"
    }
    data2plot <- data2plot[lib == sName, lapply(.SD, sum), by = VpJ, .SDcols = "count"][,interval := unlist(lapply(count, f))]
    breaks <- unique(data2plot[, interval])
    plotBreaks <- breaks[order(nchar(breaks), breaks)]
    data2plot <- data2plot[, lapply(.SD, sum), by = interval, .SDcols = "count"][,percent := prop.table(count)]
    p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = interval, y = percent, fill = interval)) +
            ggplot2::geom_bar(stat = "identity", show.legend=F) +
            ggplot2::scale_x_discrete(limits=plotBreaks) +
            ggplot2::ylim(0, 1) +
            ggplot2::geom_text(ggplot2::aes(y = percent, label = paste(100*round(percent, 2), "%")), vjust=-1, size=7) +
            ggplot2::labs(title = "TR sequence distribution", x = "count", y = "# of occurences") +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), 
                text = ggplot2::element_text(size=18), 
                axis.text.x = ggplot2::element_text(size=16),
                axis.text.y = ggplot2::element_text(size=16))
    print(p)
}

#' plot Spectratype by V genes
#'
#' function produces barplot according to the CDR3 length for each V gene.
#' @param x an object of class RepSeqExperiment.
#' @param sampleName a string indicating the sample to plot. sampleName must be existing in the slot "sampleData".
#' @param scale the type of bars in term of counts, percentages or count per million.
#' @param showCDR3 blabla
#' @return a wrapped of barplots
#' @export
# @example
plotSpectratypeV <- function(x, sampleName = NULL, scale = c("counts", "percent", "cpm"), showCDR3 = F) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
    cpm=percent=N=CDR3length=pep=lib=CDR3aa <- NULL
    sNames <- rownames(sData(x))
    if (is.null(sampleName)) {
        index <- sNames[1]
        #cat("Plot for the first sample in x:", index, ".\n")
    } else {
        index <- sampleName 
        if (length(sampleName) > 1) {
            index <- sampleName[1]
            #cat("Only the first sample is used, sample:", index, ".\n")           
        }                  
        if (is.na(match(index, sNames))) {
            stop("Sample ", index, " not found in x.")                 
        }
    }
    scl <- match.arg(scale)
    cts <- data.table::copy(assay(x))
    # subset count
    if (scl == "counts"){
        tmp <- cts[lib == index,][, CDR3length := nchar(CDR3aa)]
        data2plot <- tmp[, .(.N), by = .(V, CDR3length)]
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = N))
    }
    if (scl == "percent"){
        data2plot <- cts[lib == index,][, CDR3length := nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,percent := prop.table(N), by = V]#geom_bar
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = percent)) + 
            ggplot2::scale_y_continuous(labels = scales::percent)
    }
    if (scl == "cpm"){
        data2plot <- cts[lib == index,][, CDR3length:=nchar(CDR3aa)][,.(.N), by = .(V, CDR3length)][,cpm := prop.table(N)*10^6, by = V]#geom_bar
        p <- ggplot2::ggplot(data = data2plot, ggplot2::aes(x = CDR3length, y = cpm))
    }
    p <- p + ggplot2::geom_bar(stat = "identity") + 
            ggplot2::scale_x_continuous(breaks = data2plot[, unique(CDR3length)]) +
            ggplot2::labs(title = paste(index,": V distribution by aa length"), x = "CDR3 length (aa)", y = scl) +
            ggplot2::facet_wrap(~ factor(V, levels = naturalsort::naturalsort(unique(V))), ncol = 4, scales = "free") + 
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14), 
                text = ggplot2::element_text(size=14), 
                axis.text = ggplot2::element_text(size=14),
                axis.text.x = ggplot2::element_text(size=14),
                axis.text.y = ggplot2::element_text(size=14))
  return(p)
}


# plotSimilarityMatrix <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa"), colorBy = NULL) {
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#   variable <- NULL
#   levelChoice <- match.arg(level)
#   cols <- c("lib", levelChoice, "count")
#   tmp <- assay(x)[,..cols]
#   sdata <- sData(x)
#   if(!is.null(colorBy)){
#   if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
#   sNames <- unique(sdata[, colorBy])
#   tmp[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#   data <- dcast(data = tmp, paste(levelChoice, "~", colorBy), value.var = "count")
#   matrixData <- as.matrix(data[,..sNames])
#   simmat <- proxy::simil(x = matrixData, by_rows = F, diag = T, upper = T)
#   mat <- as.matrix(simmat, diag = 0)
#   diag(mat) <-1
#   graph.title <- paste0(colorBy," similarity heatmap : ", levelChoice)
#   if (requireNamespace("pheatmap", quietly = TRUE)) {
#     p = pheatmap::pheatmap(mat, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 15)#, cellwidth = 12, cellheight = 12)
#   }
#   }
#   else{
#     sNames <- rownames(sdata)
#     data <- dcast(data = tmp, paste(levelChoice, "~lib"), value.var = "count")
#     matrixData <- as.matrix(data[,..sNames])
#     simmat <- proxy::simil(x = matrixData, by_rows = F, diag = T, upper = T)
#     mat <- as.matrix(simmat, diag = 0)
#     diag(mat) <-1
#     graph.title <- paste(" similarity heatmap : ", levelChoice)
#     if (requireNamespace("pheatmap", quietly = TRUE)) {
#       p = pheatmap::pheatmap(mat, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 15)#, cellwidth = 12, cellheight = 12)
#     }
#   }
#   return(p)
# }

#' create heatmap of pairwise distance matrix between samples 
#'
#' function plot a heatmap of a squared distance matrix computed between samples
#' @param x an object of class RepSeqExperiment
#' @param level a string of characters indicating the level at which performed all calculations.
#' @param method distance computation method. Choose one among manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, mahalanobis.
#' @param binary if TRUE data will be transformed into present/absent data.
#' @return a heatmap
#' @export
# @example
plotDissimilarityMatrix <- function(x, level = c("VpJ", "V", "J", "VJ", "CDR3aa"), 
                                    method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", 
                                    "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", 
                                    "binomial", "chao", "cao", "mahalanobis"), 
                                    binary = "FALSE") {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    variable <- NULL
    levelChoice <- match.arg(level)
    methodChoice <- match.arg(method)
    cols <- c("lib", levelChoice, "count")
    tmp <- data.table::copy(assay(x))[, ..cols]
    sdata <- sData(x)
    sNames <- rownames(sdata)
    groups <- sdata[, unlist(lapply(sdata, is.factor)), drop = F]
    dat <- data.table::dcast(data = tmp, paste(levelChoice, "~lib"), value.var = "count", fun.aggregate = sum)
    simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag = TRUE, upper = TRUE, binary = as.logical(binary)), .SDcols=sNames]
    if (requireNamespace("pheatmap", method = methodChoice, quietly = TRUE)) {
        p <- pheatmap::pheatmap(simmat, main = paste0("dissimilarity heatmap : ", levelChoice), 
            cluster_rows = TRUE, cluster_cols = TRUE,
            treeheight_row = 0L, clustering_distance_rows = simmat, clustering_distance_cols = simmat, 
            annotation_col=groups, show_colnames=T, labels_col = sNames,
            show_rownames=FALSE, clustering_method = "ward.D", silent = TRUE)
    }
    return(p)
}

#' biplot Multidimensional Scaling result
#'
#' function computes pairwise distances between all libs for a chosen level. Then, multidimentional scaling (2 dimensions) is applied to the resulting matrix, ggplot was used to represente individuals on the scaling dimensions.
#' @param x an object of class RepSeqExperiment.
#' @param level repertoire level VpJ, V, J, VJ, CDR3aa.
#' @param method distance computation method. Choose one among manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao, mahalanobis.
#' @param colGrp color by group of samples. A string of character containng the name of the varibale indicating group of samples (must be a column name in the slot sampleData). If NULL, the number of colors is equal to the number of sample.  
#' @return a dot plot.
#' @export
# @example
plotMDS <- function(x, level = c("VpJ", "V", "J", "VJ", "CDR3aa"), 
                        method = c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", 
                        "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", 
                        "binomial", "chao", "cao", "mahalanobis"), colGrp = NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    variable <- NULL
    levelChoice <- match.arg(level)
    methodChoice <- match.arg(method)
    cols <- c("lib", levelChoice, "count")
    tmp <- data.table::copy(assay(x))[, ..cols]
    sdata <- sData(x)
    if (is.null(colGrp)) {
        colGrp <- "Sample"
    } else {
        if (length(grep(colGrp, colnames(sdata))) == 0) colGrp <- "Sample" else colGrp <- colGrp
        }
    sNames <- rownames(sdata)
    groups <- sdata[,unlist(lapply(sdata, is.factor)), drop = F]
    dat <- data.table::dcast(data = tmp, paste(levelChoice, "~lib"), value.var = "count", fun.aggregate = sum)
    simmat <- dat[, vegan::vegdist(t(.SD), method = methodChoice, diag=TRUE, upper=TRUE), .SDcols = sNames]
    fit <- stats::cmdscale(simmat, k = 2)
    colnames(fit) <- c("D1" ,"D2")
    data2plot <- data.frame(merge(fit, sdata, by = 0), row.names=1)
    fact <- eval(parse(text = paste0("sdata$", colGrp)))
    ade4::s.class(fit, fac = fact, col=1:nlevels(fact), sub = "Multidimensional scaling")
}

#, sub=paste("Level:", levelChoice, "; Distance:", methodChoice)
# plotFrequencySpectrum <- function(x, groupBy = FALSE, colorBy = NULL){
#   if (missing(x)) stop("x is missing.")
#   if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#   if(!is.null(colorBy)){
#     if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
#     sdata <- sData(x)
#     counts <- assay(x)
#     if(groupBy){
#       counts[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#       counts <- counts[,lapply(.SD, sum), by = c("VpJ", paste(colorBy)), .SDcols = "count"][,.N, by = c(paste(colorBy), "count")]
#       print(head(counts))
#       ggplot(data = counts, aes_string(x = "count", y = "N", colour = paste(colorBy))) + geom_point() +geom_line()+ scale_x_log10()
#     }else{
#       counts <- counts[,lapply(.SD, sum), by = c("VpJ", "lib"), .SDcols = "count"]
#       counts[,paste(colorBy) := lapply(.SD, function(x){sdata[x, colorBy]}), .SDcols = "lib"]
#       data2plot <- counts[,.N, by = c("lib", "count", paste(colorBy))]
#       return(ggplot(data = data2plot, aes_string(x = "count", y = "N", colour = paste(colorBy))) + geom_point() +geom_line(aes(group = lib))+ scale_x_log10())
# 
#     }
#   }
#   else{
#     counts <- assay(x)[,lapply(.SD, sum), by = .(VpJ, lib), .SDcols = "count"][,.N, by = .(lib, count)]
#     ggplot(data = counts, aes(x = count, y = N, colour = lib)) + geom_point() +geom_line()+ scale_x_log10()
#   }
# }

#' plot of frequency 
#' 
#' function plots frequency of 
#' @param x x
#' @param colorBy color by factor name
#' @param groupBy boolean if summary by factor or not
#' @return a graph
#' @export
# @example
plotFrequencySpectrum <- function(x, colorBy = NULL, groupBy = FALSE) {
    N=std=group <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    counts <- cts[, lapply(.SD, sum), by = .(VpJ, lib), .SDcols = "count"][,.N, by = .(lib, count)]
    if (!is.null(colorBy)) {
        if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
        sdata <- sData(x)
        counts[, group := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "lib"]
        if (groupBy) {
            counts[, std := lapply(.SD, sd), by = c("count", "group"), .SDcols = "N"]
            counts <- counts[, lapply(.SD, mean), by = c("group", "std", "count"), .SDcols = "N"]
            #print(counts)
            p <- ggplot2::ggplot(data = counts, ggplot2::aes(x = count, y = N, colour = group)) + 
                    ggplot2::geom_errorbar(ggplot2::aes(ymin = N-std, ymax = N+std), width=.1, position = "dodge") + 
                    ggplot2::geom_point() +
                    ggplot2::geom_line() + 
                    ggplot2::scale_x_log10()
        } else {
            #print(counts)
            p <- ggplot2::ggplot(data = counts, ggplot2::aes(x = count, y = N, colour = group)) + 
                    ggplot2::geom_point() +
                    ggplot2::geom_line(ggplot2::aes(group = lib)) + 
                    ggplot2::scale_x_log10()
        }
    } else  {
        p <- ggplot2::ggplot(data = counts, ggplot2::aes(x = count, y = N, colour = lib)) + 
                ggplot2::geom_point() +
                ggplot2::geom_line() + 
                ggplot2::scale_x_log10()
    }
    return(p)
}

# plot rarefaction curve
#
# function plot rarefaction curves
# @param x
# @param colorBy
# @param groupBy
# @return a graph 
#plotRarefaction <- function(x, colorBy = NULL, groupBy = FALSE){
#    if (missing(x)) stop("x is missing.")
#    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
#    sdata <- sData(x)
#    #print(colorBy)
#    cts <- data.table::copy(assay(x))
#    if (is.null(colorBy)) {
#        dcounts <- data.table::dcast(cts, VpJ~lib, value.var = "count", fun.aggregate = sum)
#        #cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[, !"VpJ"])]
#    } else {
#        sNames <- unique(sdata[, colorBy])
#        if (is.na(match(eval(colorBy), colnames(sdata)))) stop(paste0(colorBy," not found in sData(x)."))
#        if (groupBy) {
#            counts <- cts[, paste(colorBy) := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "lib"]
#            dcounts <- data.table::dcast(data = counts, paste("VpJ~",colorBy), value.var = "count", fun.aggregate = sum)
#            #print(dcounts[,!"VpJ"])
#            #cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[,!"VpJ"])]
#        } else {
#            counts <- assay(x)
#            dcounts <- dcast(data = counts, VpJ~lib, value.var = "count", fun.aggregate = sum)
#            #cumulcounts <- dcounts[, lapply(.SD, cumsum), .SDcols = colnames(dcounts[,!"VpJ"])]
#        }
#    }
#    #print(nbcols(counts[, lapply(.SD, function(x){x==0}), .SDcols =  colnames(dcounts[,!"VpJ"])]))
#    # VpJlist <- as.vector(counts[, VpJ])
#    #print(dcounts)
#    rare <- rarefaction(dcounts[, !"VpJ"], 10)#, step)
#    nblines <- ncol(rare)
#    raredt <- data.table::setDT(rare)
#    #raredt[, size := c(1, step*1:(nrow(raredt)-1))]
#    rarem <- melt(raredt, measure.vars = colnames(raredt[,!"size"]))
#    if (!is.null(colorBy) & !groupBy) {
#        rarem[, paste(colorBy) := lapply(.SD, function(x) sdata[x, colorBy] ), .SDcols = "variable"]
#        p <- ggplot2::ggplot(data = rarem[value!=0], aes_string(x = "size", y = "value", colour = colorBy)) + 
#                ggplot2::geom_line()
#        return(p)
#    }
#    ggplot(data = rarem[value!=0], aes(x = size, y = value, colour = variable)) + geom_line(size = 2)
#  # tcounts = t(as.matrix(dcounts[,!"VpJ"]))
#  # #print(head(tcounts))
#  # raremax <- min(rowSums(tcounts))
#  # rarecurve(tcounts, step = 10000, sample = raremax,
#  #            label = FALSE)
#}


#' function plotDistribVpJ
#'
#' function 
#' @param x an object of class RepSeqExperiment
#' @param colorBy column name in sample Data indicating groups of samples
#' @param aggreg aggregation of counts by sum or by mean
#' @return a plot
#' @export
# @example
plotDistribVpJ <- function(x, colorBy = NULL, aggreg = c("sum", "mean")){
    group <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if(is.null(colorBy)) stop("need group for now")
    sdata <- sData(x)
    counts <- data.table::copy(assay(x))
    counts[, group := lapply(.SD, function(x) sdata[x, colorBy] ), .SDcols = "lib"]
    aggr <- match.arg(aggreg)
    if (aggr == "sum") counts <- counts[, lapply(.SD, sum), by = .(group, VpJ), .SDcols = "count"]
    if (aggr == "mean") counts <- counts[, lapply(.SD, mean), by = .(group, VpJ), .SDcols = "count"]
    counts[, rank := lapply(.SD, frankv, ties.method = "min", order = -1L), by = group, .SDcols = "count"]
    counts <- unique(counts[,!"VpJ"])
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = rank, y = count, colour = group)) + 
           ggplot2::geom_point() + 
           ggplot2::scale_x_log10() + 
           ggplot2::scale_y_log10() + 
           ggplot2::ylab(paste(aggr, "counts")) +
           ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14), 
                text = ggplot2::element_text(size=14), 
                axis.text = ggplot2::element_text(size=14),
                axis.text.x = ggplot2::element_text(size=14),
                axis.text.y = ggplot2::element_text(size=14))
    return(p)
}

#' plot Venn's diagram
#'
#' function plot Venn's diagram up to 5 samples
#' @param x an object of class RepSeqExperiment
#' @param level repertoire level
#' @param colorBy color by group of samples. If not NULL 
#' @return a Venn's diagram
#' @export
# @example
plotVenn <- function(x, level = c("V", "J", "VJ", "VpJ", "CDR3aa"), colorBy = NULL) {
    group <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    sdata <- sData(x)
    counts <- data.table::copy(assay(x))
    if (length(colorBy) > 1) {
        sNames <- colorBy
        dcounts <- data.table::dcast(data = counts, paste(levelChoice,"~lib"), value.var = "count", fun.aggregate = sum)
    } else {
        sNames <- unique(sdata[, colorBy])
        counts[, group := lapply(.SD, function(x) sdata[x, colorBy]), .SDcols = "lib"]
        dcounts <- data.table::dcast(data = counts, paste(levelChoice,"~group"), value.var = "count", fun.aggregate = sum)
    }
    a <- limma::vennCounts(dcounts[ , ..sNames])
    limma::vennDiagram(a, circle.col = c("red", "blue", "green3", "purple", "yellow"))
}

#' plot heatmap of counts or percentages 
#'
#' function plots a heatmap of counts or percentages of sorted 
#' @param x an object of class RepSeqExperiment
#' @param level of the repertoire
#' @param type count or usage
#' @return a heatmap 
#' @export
# @example
plotmuScore <- function(x, level=c("V", "J", "VJ"), type=c("count", "usage")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    typeChoice <- match.arg(type)
    sdata <- sData(x)
    #grp <- which(unlist(lapply(sdata, nlevels)) < nrow(sdata))
    #sdata <- sdata[, grp, drop = FALSE]
    # replace of forbidden characters
    #sdata <- apply(sdata, 2, function(tab) gsub("\ ", ".", tab)) # replace blank characters by "dot"
    #sdata <- apply(coldata, 2, function(tab) gsub("\\+", "p", tab)) # replace + by "p"
    #sdata <- apply(coldata, 2, function(tab) gsub("\\-", "m", tab)) # replace - by "m" 
    #sdata <- data.frame(coldata) # replace - by "m"
    #rownames(coldata) <- gsub("\\-", ".", rownames(sdata))
    sNames <- rownames(sdata)
    groups <- sdata[, unlist(lapply(sdata, is.factor)), drop = F]
    data2plot <- data.frame(muScore(x, levelChoice, typeChoice), row.names=1)
    p <- pheatmap::pheatmap(as.matrix(data2plot[, rownames(groups), drop=FALSE]), 
            main = paste(typeChoice, "of", levelChoice, "in samples"), cluster_rows = FALSE, cluster_cols = TRUE,
            treeheight_row = 0L, annotation_col=groups, show_colnames=T, labels_col = sNames,
            show_rownames=TRUE, clustering_method = "ward.D2", silent = TRUE, color = grDevices::colorRampPalette(c("lightgrey", "red"))(100))
    return(p)
}

#' sample against sample plot
#'
#' function plots 
#' @param x an object of class RepSeqExperiment
#' @param level choose between V, J, VJ, VpJ, CDR3aa
#' @param libs a vector of characters of length 2 containing the sample names to plot.
#' @param scale normal or log10 scale of axes
#' @return a dot plot
#' @export
# @example
plot2v2count <- function(x, level = c("V", "J", "VJ", "VpJ", "CDR3aa"), libs = NULL, scale = c("counts", "log")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    scaleChoice <- match.arg(scale)
    cols <- c("lib", levelChoice, "count")
    cts <- data.table::copy(assay(x))
    counts <- cts[lib ==libs[1] | lib == libs[2], ..cols]
    data2plot <- data.table::dcast(counts, paste(levelChoice, "~lib"), value.var="count", fun.aggregate = sum)
    p <- ggplot2::ggplot(data2plot, ggplot2::aes_string(libs[1], libs[2])) +
        ggplot2::geom_count(ggplot2::aes(alpha = ..n..), size = 2) + 
        ggplot2::scale_alpha(range = c(0.25, 1)) +
        ggplot2::labs(title = paste("counts comparison per ", levelChoice)) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
    p2 <- p + ggplot2::scale_x_log10() + 
                ggplot2::scale_y_log10() + 
                ggplot2::annotation_logticks()
    if (scale == "log") return(p2) else return(p)
}


