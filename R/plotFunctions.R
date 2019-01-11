#------------------- plot functions --------------------#
#' plot a heatmap og V-J count 
#'
#' function plots a heatmap of V-J counts for sampleName 
#'
#' @param x a matrix of count data
#' @param sampleName sample to plot, sample name must be existed among x column names.
#' @return a heatmap of count
#' @export
#' @seealso \code{\link[pheatmap]{pheatmap}}
# @example
# plotCountVJ()
plotCountVJ <- function(x, sampleName=NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
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
    tmp <- dcast(assay(x)[lib==index, ], V~J, value.var="count", fun.aggregate = sum)
    data2plot <- data.frame(tmp, row.names=1)
    graph.title <- paste0("sample: ", index)
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(data2plot, main=graph.title, cluster_rows = FALSE, cluster_cols = FALSE, heigh = 3)
    }
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
#' function plots spectratype for V genes of a repertoire.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param sampleName sample to plot, if NULL the first sample of x is plotted.
#' @return a barplot
#' @export
# @example
# plotSpectratype
plotSpectratype <- function(x, sampleName=NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is epxected.") 
    CDR3length=pep=lib=CDR3aa <- NULL
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
    # subset count
    tmp <- assay(x)[lib == index, ][, CDR3length:=nchar(CDR3aa)]
    #tmp[tmp[, .I[.SD>0], .SDcols=eval(index)]]
    # V distribution by aa length 
    pastek <- table(tmp$V, tmp$CDR3length)
    # compute percentages
    #pastek <- prop.table(pastek, margin=1)
    # plot
    par(xpd=TRUE, mar=c(5.1, 4.1, 4.1, 9.1))
    barplot(pastek, main=paste(index,": V Distribution by aa length"), xlab="CDR3 length (aa)", col=1:nrow(pastek)) 
    legend("topright", inset=c(-0.2, 0), rownames(pastek), col=1:nrow(pastek), pch=19, cex=0.5)  
    par(xpd=FALSE, mar=c(5.1, 4.1, 4.1, 2.1))
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
plotRenyiProfiles <- function(x, alpha=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=c("VpJ", "V", "J", "VJ", "CDR3aa"), colorBy=NULL) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (length(alpha) < 2) stop("At least 2 values of alpha is needed.")
    variable <- NULL
    sNames <- rownames(sData(x)) 
    levelChoice <- match.arg(level)
    # compute Renyi
    data2plot <- renyiProfiles(x, scales=alpha, level=levelChoice)
    # plot 
    if (is.null(colorBy)) {
        matplot(as.matrix(data2plot[, variable]), as.matrix(data2plot[, ..sNames]), type="l", 
            xlab="alpha", ylab="Renyi's entropy", main=paste0("Level ", levelChoice,": Renyi profiles"))
        } else {
            if (is.na(match(eval(colorBy), colnames(sData(x))))) stop(paste0(colorBy," not found in sData(x)."))
            # get factor
            fact <- sData(x)[, eval(colorBy)]
            if (!is.factor(fact)) stop(paste0(colorBy, "musts be a factor."))
            # convert factor to color
            color <- as.numeric(fact)    
            # plot
            par(xpd=TRUE, mar=c(5.1, 4.1, 4.1, 8.1))
            matplot(as.matrix(data2plot[, variable]), as.matrix(data2plot[, ..sNames]), type="l", 
             xlab="alpha", ylab="Renyi's entropy", main=paste("Level ", levelChoice,": Renyi profiles"), col=color)
            legend("topright", inset=c(-0.2, 0), levels(fact), col=unique(color), pch=19, cex=0.7)  
            par(xpd=FALSE, mar=c(5.1, 4.1, 4.1, 2.1))
        }
}









