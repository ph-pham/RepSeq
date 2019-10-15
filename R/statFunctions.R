#---------------- statistics -------------#

#' RepSeqCount to DESeqDataSet 
#' 
#' function transforms a RepSeqCount object to an DESeqDataSet object
#'
#' @param x an object of class RepSeqExperiment
#' @param conditions name(s) of column(s) in sample information (sData(x)) defining test between biological conditions.
#' @param level of repertoire to analyze \code{VpJ, VJ, V, J, CDR3aa}.
#' @return a DESeqDataSet object
#' @export
# @example
# example here
toDESeq2 <- function(x, conditions, level=c("VpJ", "V", "J", "VJ", "CDR3aa")) {
    # test 
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("a RepSeqCount object is expected.")
    if (missing(conditions)) stop("user musts provide at least 1 condition.")    
    # get sample info
    coldat <- sData(x)[, conditions, drop=FALSE]
    coldat <- apply(coldat, 2, function(x) gsub("\ ", ".", x)) # replace blank characters by "dot"
    coldat <- apply(coldat, 2, function(x) gsub("\\+", "p", x)) # replace + by "p"
    coldat <- apply(coldat, 2, function(x) gsub("\\-", "m", x)) # replace - by "m"
    levelChoice <- match.arg(level)
    cts <- countFeatures(x, level=levelChoice)
    if (length(conditions) > 1) conditions <- paste(conditions, collapse="+")
    rownames(coldat) <- gsub("-", ".", rownames(coldat))
    cts <- data.frame(cts, row.names=1)
    colnames(cts) <- gsub("-", ".", colnames(cts))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = coldat, design = as.formula(paste0("~" ,conditions)))
    return(dds)
}


#' estimate size factor`
#'
#' function 
#'
#' @param x an object of class RepSeqExperiment
#' @param level level of repertoire to analyze \code{VpJ, VJ, V, J, CDR3aa}.
#' @param method normaliztion method used for size factor computation.
#' @param UsePseudoRef a boolean indicating if Chao indices will be computed according to a reference repertoire (geometric mean repertoire across all samples). 
#' @return a vector of normalized size factors
#' @export
# @example
estimateSF <- function(x, level=c("VpJ", "CDR3aa"), method=c("Chao", "iChao", "worChao", "Chao.gmmean", "Chao.median"), UsePseudoRef=TRUE) {
    # test 
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("a RepSeqCount object is expected.")
    chaoest <- function(x) {
        return(.chao1(x)$chao.est)
    }
    gm_mean <- function(x, na.rm=TRUE) {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    levelChoice <- match.arg(level)
    meth <- match.arg(method)
    dat <- copy(assay(x))
    chao1 <- ichao <- chaowor <- sj <- s0 <- NULL
#    if (meth == "GMPR") {
#        cts <- countFeatures(x, level=levelChoice)
#        out <- GMPR::GMPR(t(data.frame(cts, row.names=1)))
#        }
    if (meth == "Chao") {
        chao <- dat[, .(count=sum(count)), by=c("lib", levelChoice)][, .(s0=sum(count>0), sj=sum(count), chao1=chaoest(count)), by="lib"]
        if (UsePseudoRef) {
            ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
            chao.ref <- chaoest(ref$gmmean)
            out <- chao[, (chao.ref/s0) * sj]
            } else out <- chao[, (chao1/s0) * sj]
        names(out) <- chao$lib
        }
    if (meth == "Chao.gmmean") {
        chao <- dat[, .(count=sum(count)), by=c("lib", levelChoice)][, .(chao1=chaoest(count)), by="lib"]
        out <- dat[, .(s0=sum(count>0), sj=sum(count)), by="lib"][, (gm_mean(chao$chao1)/s0)*sj]
        names(out) <- tmp$lib
        }
    if (meth == "Chao.median") {
        chao <- dat[, .(count=sum(count)), by=c("lib", levelChoice)][, .(chao1=chaoest(count)), by="lib"]
        tmp <- dat[, .(s0=sum(count>0), sj=sum(count)), by="lib"][, stats::median(chao$chao1[chao$chao1>0], na.rm=TRUE)*sj/s0, by="lib"]
        out <- tmp$V1
        names(out) <- tmp$lib
        }
    if (method == "iChao") {
        chao <- dat[, .(count=sum(count)), by=c("lib", levelChoice)][, .(s0=sum(count>0), sj=sum(count), ichao=iChao(count)), by="lib"]
        if (UsePseudoRef) {
            ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
            chao.ref <- iChao(ref$gmmean)
            out <- chao[, (chao.ref/s0) * sj]
            } else out <- chao[, (ichao/s0) * sj]
        names(out) <- chao$lib
        }
    if (method == "worChao") {

        chao <- dat[, .(count=sum(count)), by=c("lib", levelChoice)][, .(s0=sum(count>0), sj=sum(count), chaowor=Chaowor(count)), by="lib"]
        if (UsePseudoRef) {
            ref <- dat[, .(gmmean=gm_mean(count)), by=levelChoice]
            chao.ref <- Chaowor(ref$gmmean)
            out <- chao[, (chao.ref/s0) * sj]
            } else out <- chao[, (chaowor/s0) * sj]
        names(out) <- chao$lib
        }
    return(out)
}


#' compute multivariate score
#'
#' function computes multivariate score for each feature at level of the repertoire \code{VJ, V, J}. The level VpJ is not supported for the moment.
#'
#' @param x a RepSeqExperiment object
#' @param level of repertoire to analyze \code{VJ, V, J, VpJ, CDR3aa}.
#' @param type type of data to compute muScore, \code{count} or \code{usage}.
#' @return a data.table ordered by the multivariate decreasing scores. The last column contains the multivariate score.
#' @export
# @example
# example here
muScore <- function(x, level=c("V", "J", "VJ", "VpJ", "CDR3aa"), type=c("count", "usage")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("a RepSeqExperiment object is expected.")
    score <- NULL
    # match arguments
    levelChoice <- match.arg(level)
    typeChoice <- match.arg(type)
    # get sample names
    sNames <- rownames(sData(x))
    if (levelChoice %in% c("VpJ", "CDR3aa")) x <- getOverlaps(x)
    # select type of data
    if (typeChoice == "usage") {
        DT <- segmentUsage(x, levelChoice)
    } else {
        DT <- countFeatures(x, levelChoice)
        }
    DT[, score := muScoreOB(.SD), .SDcols=sNames]
    out <- DT[order(score, decreasing=TRUE)]
    return(out)
}

# old method        
    # long
#    DT.melt <- melt(DT, id.var=levelChoice, fill=0)
    # create row & column Ids
#    DT.melt[, c("colId", "rowId") := list(rep(1:length(sNames), each = nrow(DT)), rep(1:nrow(DT), times = length(sNames)))]
    
    # key for organizing
#    setkey(DT.melt, rowId)

    # computes pairwise comparison between rows (for each column, return -1 if row x < row y; 0 of not comparable; 1 if row x > row y) 
#    diff.DT <- DT.melt[
#        , {
#            ccl <- unique(rowId)
#            vv <- value
#            DT.melt[rowId != ccl, .(row2 = rowId, 1 * (vv > value) - 1 * (value > vv))]
#            }
#        , keyby = .(row1 = rowId)
#    ]
    # set keys
#    setkey(diff.DT, row1, row2)
    # sum by 
#    sum.DT <- diff.DT[, sign(sum(V2)), by=c("row1", "row2")]
    # compute score
#    score <- sum.DT[, sum(V1), by="row1"]
    # update
#    DT[, score:=score$V1]
#    out <- DT[order(score, decreasing=TRUE)]
#    return(out)
#}

# multivariate rank test
#
# function tests  
#
# @param x a data.frame or data.table contains parameter outcomes (rows) measured for samples (columns).
# @param y a data.frame or data.table contains parameter outcomes (rows) measured for samples (columns).
# @param alternative 
# @param
# @export
# @example

#mimi.test <- function(x, y, alternative="both", method=c("obrien", "wittkowski", "comb")) {
#    if (missing(x)) stop("x is missing.")
#    if (missing(y)) stop("y is missing.")
#    if (nrow(x) != nrow(y)) stop("number of parameters must agree (number of rows).")
#    n <- ncol(x)
#    m <- ncol(y)
#    dat <- setDT(cbind(x, y))
    # create list of combinations
#    comblist <- combn(c(colnames(x), colnames(x)), 2, FUN=list)
    # compute differences
#    diff.cols <- dat[, lapply(comblist, function(x) rij(get(x[1]), get(x[2]))) ]
    # 
    #setnames(diff.cols, names(diff.cols), sapply(comblist, paste, collapse="_"))
#    u.method <- match.arg(method)
#    if (u.method =="obrien") {
        # compute U score
#        U <- diff.cols[, lapply(.SD, sum)][, Reduce(`+`, .SD)]
#        U <- U/(n * m)
#        }
#    if (u.method == "wittkowski") {
#        U <- diff.cols[, lapply(.SD, function(x) 1*(sum(x)>0) - 1*(sum(x)<0))][, Reduce(`+`, .SD)]
#        U <- U/(n * m)
#    }
    
#return(res)
#}

#' create contrast matrix from factor level names
#'
#' function creates contrast matrix from factor level names. 
#'
#' @param x a factor
#' @param ... other parameters will be passed to function \code{contr.sum}
#' @return a contrast design matrix
#' @export 
named.contr.sum <- function(x, ...) {
    if (is.factor(x)) {
        x <- levels(x)
    } else if (is.numeric(x) & length(x)==1L) {
        stop("cannot create names with integer value. Pass factor levels")
    }
    x <- stats::contr.sum(x, ...)
    colnames(x) <- apply(x,2,function(x) 
         paste(names(x[x>0]), names(x[x<0]), sep="-")
    )
    x
}

#' Get normalized count 
#'
#' function computes the estimated size factor according to the choice and level of the repertoire, and the  
#'
#' @param x a RepSeqExperiment. 
#' @param method method used for normalization. 
#' @param UsePseudoRef a boolen indicatif whether a reference repertoire will be used for normalizaition. 
#' @return an object of class RepSeqExperiment with normalized counts.
#' @export
# @example
# Discussion https://support.bioconductor.org/p/66067/
normalizeCounts <- function(x, method=c("Chao", "iChao", "worChao", "Chao.gmmean", "Chao.median"), UsePseudoRef=TRUE) {
    # test 
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("a RepSeqCount object is expected.")
    choice <- match.arg(method)
    dat <- copy(assay(x))
    sampleinfo <- sData(x) 
    sf <- estimateSF(x, level="VpJ", method=choice, UsePseudoRef=UsePseudoRef)
    dat[, count:=count/rep(sf, table(dat$lib))]
    sampleinfo$sf <- sf
    x.hist <- data.frame(rbind(History(x), history = paste0("normalizedCounts; x=", deparse(substitute(x)), "; method=", choice, "; UsePseudoRef=", UsePseudoRef)), stringsAsFactors = FALSE)
    out <- methods::new("RepSeqExperiment", assayData=dat, sampleData=sampleinfo, metaData=mData(x), History=x.hist) 
}




