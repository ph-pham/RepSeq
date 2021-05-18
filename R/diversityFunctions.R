#-------- diversity computation -----------#

#' updates sample data with basic statistics
#'
#' function produces several useful diversity indices. 
#'
#' @param x an object of class RepSeqExperiment
#' @param level level of reppertoire to assess diversity, choice are VpJ, V, J, VJ or CDR3aa)
#' @return a data frame of Chao
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' divindices <- basicIndices(RepSeqData, level="VpJ")
#' divindices
#' }
basicIndices <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa")) {
    if (missing(x)) stop("x is missing. An object of class RepSeqExperiment is expected.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    # summary function
    pastek <- function(y) list(shannon=.diversity(y, index="shannon"), simpson=.diversity(y, index="simpson"), 
                            invsimpson=.diversity(y, index="invsimpson"), gini=.gini(y), 
                            chao1=.chao1(y)$chao.est, chao1.se=.chao1(y)$chao.se, ichao=iChao(y), chaowor=Chaowor(y))   
    out <- copy(assay(x))[, .(count=sum(count)), by=c("lib", levelChoice)][, pastek(count), by="lib"]
    return(out)
}

# function used internally
.diversity <- function(x, index=c("shannon", "simpson", "invsimpson"), norm=FALSE, base=exp(1)) {
    if (missing(x)) stop("data set x is required.")
    x <- x/sum(x)
    id <- match.arg(index)
    if (id == "shannon") {
        x <- -x * log(x, base)
        } else {
            x <- x * x
            }
    H <- sum(x, na.rm = TRUE)
    if (norm) H <- H/log(sum(x>0), base)
    if (id == "simpson") {
        H <- 1 - H
        } else if (id == "invsimpson") {
            H <- 1/H
            }
    return(H)
}

#' compute diversity indices
#'
#' function computes diversity indices for level \code{VJ}, \code{V} or \code{J}
#'
#' @param x a object of class RepSeqExperiment.
#' @param index diversity index, the choice are "shannon", "simpson", "invsimpson".
#' @param level diversity index should be computed for the level "V", "J", or "VJ".
#' @param norm normalized diversity index will be returned (divided by log(S) with S number of observed species).
#' @param base log base . Default value: exp(1)
#' @return a data.table of diversity indices computed at the level "level".
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' divindices <- divLevel(x = RepSeqData, index = "shannon", level = "VJ")
#' divindices
#' }
divLevel <- function(x, index = c("shannon", "simpson", "invsimpson"), level = c("VJ", "V", "J"), norm = TRUE, base = exp(1)) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
    divIndex <- match.arg(index)
    levelChoice <- match.arg(level)
    res <- assay(x)[, .diversity(count, index=index, base=base, norm=norm), by=c("lib", levelChoice)]
    out <- data.table::dcast(res, as.formula(paste0(levelChoice, "~lib")), value.var="V1", fill=0)
    return(out)
}

# compute Renyi or hill number
#
# function computes Renyi (Hill's numbers) indices according to \alpha.
# @param x a vector of counts.
# @param alpha a number between 0 and infinity, alpha is Renyi's parameter.
# @param hill a boolean if TRUE the Hill's indice is computed.
# @return a number
# @export
.renyiCal <- function(x, alpha=2, hill = FALSE) {
    if (missing(x)) stop("x is missing.")
    if (!is.numeric(x)) stop("x musts be a numerical vector.")
    x <- x/sum(x)
    if (alpha != 0 && alpha != 1 && alpha != Inf) {
        res <- log(sum(x^alpha))/(1 - alpha)
        } else {
            if (alpha == 0) {
                res <- log(sum(x > 0))
            } else if (alpha == Inf) {
                    res <- -log(max(x))
                    } else {
                        res <- sum(-x * log(x, exp(1)), na.rm=TRUE)
                        }
            }
    if (hill) res <- exp(res)
    return(res)
}

#' compute Renyi's profile
#'
#' function computes Renyi profile for all samples in a RepSeqExperiment object
#'
#' @param x an object of class RepSeqExperiment
#' @param scales a numerical vector of alpha to evaluate Renyi's entropy
#' @param level level of repertoire to assess Renyi profiles
#' @return a data.table containing Renyi's profiles computed for each repertoire 
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' renyiNumbers <- renyiProfiles(RepSeqData, level = "VpJ")
#' renyiNumbers
#' }
renyiProfiles <- function(x, scales = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level = c("VpJ", "V", "J", "VJ", "CDR3aa")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is exprected.")
    levelChoice <- match.arg(level)
    out <- copy(assay(x))[, .(count=sum(count)), by=c("lib", levelChoice)][, lapply(scales, function(y) .renyiCal(count, y)), by="lib"]
    data.table::setnames(out, c("lib", scales))
    out <- data.table::dcast(melt(out, id.vars = "lib"), variable ~ lib)
    return(out)
}

# compute Gini's coefficient
#
#
#
.gini <- function(x) {
    x <- sort(x)
    n <- length(x)
    out <- 1/n * (n + 1 - 2 * sum((n + 1 - 1:n) * x)/sum(x))
    out
}

# compute Chao1 indices
#
# function computes Chao indices
#
# @param x a vector of counts
# @return a list containing estimated Chao1's diversity and its standard error.
.chao1 <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    S <- sum(x > 0)
    #if (f2 > 0) est <- S + f1^2/(2*f2)
    #if (f2 == 0) est <- S + f1*(f1 - 1)/2
    est <- S + f1*(f1 - 1)/(2 * (f2 + 1))
    r <- f1/f2
    chao.var <- f2 * ((r^4)/4 + r^3 + (r^2)/2)
    chao.se <- sqrt(chao.var)
    #return(est)
    return(list(chao.est=est, chao.se=chao.se))
}

#' Improved Chao1
#'
#' function computes the improve version of Chao1
#' @param x a vector of count.
#' @return improved Chao1 value (ref) Chao and Lin Chao, A. and Lin, C.-W. (2012) Nonparametric lower bounds for species richness and shared species richness under sampling without replacement. Biometrics,68, 912–921. 
#' @export
#' @examples
#' set.seed(1234)
#' x <- rbinom(20, 5, 0.5)
#' y <- rbinom(20, 20, 0.1)
#' iChao(x)
#' iChao(y)
iChao <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    f3 <- sum(x == 3)
    f4 <- sum(x == 4)
    if (f4 == 0) f4 <- 1
    n <- sum(x)
    p1 <- (n-3)/n
    p2 <- (n-3)/(n-1)
    est <- .chao1(x)$chao.est + p1*f3/(4*f4) * max(f1 - p2*f2*f3/(2*f4), 0)
    return(est)
}

#' Adjusted Chao1 for sampling without replacement 
#'
#' function compute the correct 
#' @param x a vector of counts
#' @return a value of Chao1
#' @export
#' @examples
#' set.seed(1234)
#' x <- rbinom(20, 5, 0.5)
#' y <- rbinom(20, 20, 0.1)
#' Chaowor(x)
#' Chaowor(y)
Chaowor <- function(x) {
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    S <- sum(x > 0)
    n <- sum(x)
    q <- n/S
    dn1 <- n*2*f2/(n-1)
    dn2 <- q*f1/(q-1)
    est <- S + (f1^2)/(dn1 + dn2)
    return(est)
}

# rarefaction curve https://dave-clark.github.io/post/speeding-up-rarefaction-curves-for-microbial-community-ecology/
# compute rarefaction data
#
# function compute 
# @param x a numeric vector
# @return a rarefied vector
# @export
rarefyDT <- function(x) {
    if (missing(x)) stop("x is missing, a vector of integers is expected.")
    lib.size <- sum(x)
    breakpoint <- lib.size * 0.1
    if (breakpoint < 10) {
        steps <- pretty(1:lib.size)
        steps[length(steps)] <- lib.size
        } else {
            step1 <- pretty(1:breakpoint, n = 10)
            step2 <- pretty(breakpoint:sum(x), n = 10)
            step2[length(step2)] <- lib.size
            steps <- c(step1, step2[step2>0])
            }
    xx <- vegan::rarefy(x, sample=steps)
    output <- data.frame(x = attr(xx, "Subsample"), y=as.double(xx))
    output$x <- as.double(output$x)
    return(output)
}

#' compute rarefaction table for a RepSeqExperiment object
#'
#' function computes rarefaction table for each sample (libnames). The number of clontypes detected in the function of the number of sequences generated.
#' @param x an object of class RepSeqExperiment
#' @return a data.table containing rarefaction table in order to be represented graphically as rarefaction curves.
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' raretab <- raretabRepSeq(x = RepSeqData)
#' raretab
#' }
raretabRepSeq <- function(x) {
    if (missing(x)) stop("x is required. An object of class RepSeqExperiment is expected.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected")
    cts <- data.table::copy(assay(x))
    raretab <- cts[, rarefyDT(count), by = lib]
    return(raretab)
}

#' down-sampling to have the same library size.
#' 
#' function resamples counts to a lower same library size. Two resampling methods is proposed: with or without replacement.  
#' @param x an object of class RepSeqExperiment.
#' @param sample.size an integer indicate the final library size for all samples. The default is the smallest library size among all samples of the data set.
#' @param rngseed a integer used as seed for a reproducible result.
#' @param replace a boolean if TRUE the resampling is drawn with replacement.
#' @param verbose a boolean  if TRUE the detail of computation step will be shown.
#' @return an object of class RepSeqExperiment
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' RepSeqDataDownSampling <- rarefyRepSeqExp(x = RepSeqData, rngseed = 1234)
#' RepSeqDataDownSampling
#' }
rarefyRepSeqExp <- function(x, sample.size = min(sData(x)$nReads), rngseed = FALSE, replace = TRUE, verbose = TRUE) {
    if (missing(x)) stop("x is missing, an object of class RepSeqExperiment is expected")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected")
    sampleinfo <- sData(x)
    if (as(rngseed, "logical")) {
        set.seed(rngseed)
        if (verbose) {
            cat("`set.seed(", rngseed, ")` was applied for reproducibility of the random subsampling. Please record this number for future purpose.\n")
            cat("You can try `set.seed(", rngseed, "); .Random.seed` for the full vector.\n")
        }
    } else if (verbose) {
        cat("You set `rngseed` to FALSE. Make sure you've set & recorded\n", 
            " the random seed of your session for reproducibility.\n", 
            "See `?set.seed`\n")
    }
    if (length(sample.size) > 1) {
        warning("`sample.size` had more than one value. ", "Using only the first. \n ... \n")
        sample.size <- sample.size[1]
    }
    if (sample.size <= 0) {
        stop("sample.size = ", sample.size, " is less than or equal to zero. ", "Need positive sample size to work.")
    }
    if (sample.size > max(sampleinfo$nReads)) {
        stop("sample.size = ", sample.size, " is larger than the largest library size.")
    }
    if (verbose) {
        message("Down sampling to ", sample.size, " reads...")
    }
    if (min(sampleinfo$nReads) < sample.size) {
        lib.drop <- rownames(sampleinfo)[which(sampleinfo$nReads < sample.size)]
        if (verbose) {
            message(length(lib.drop), " samples removed", "because they contained fewer reads than `sample.size`.")
            message("Up to first five removed samples are: \n")
            message(paste(lib.drop[1:min(5, length(lib.drop))], 
                sep = "\t"))
            message("...")
        }
        x <- dropSamples(x, lib.drop)
    }
    # update assay
    cts <- data.table::copy(assay(x))
    if (replace) {
        cts[, count := {rareout = numeric(length(count)); 
                        tmp1 = sample(1:length(count), sample.size, replace = TRUE, prob = count/sum(count));
                        tmp2 = table(tmp1);
                        rareout[as(names(tmp2), "integer")] <- tmp2;
                        list(rarevec=rareout)}, by = lib]
    } else {
        cts[,count := {rareout = numeric(length(count)); 
                   tmp1 = rep(1:.N, count); 
                   tmp2 = sample(tmp1, 100, replace = FALSE); 
                   tmp3 = table(tmp2); 
                   rareout[as(names(tmp3), "integer")] <- tmp3; 
                   list(rarevec=rareout)}, by = lib]
    }
    cts <- cts[count>0]
    # update sampleData
    stats <- cts[, c(.(nReads = sum(count)), lapply(.SD, uniqueN)), by = "lib", .SDcols = c("VpJ", "V", "J", "VJ", "CDR3aa")]
    stats <- data.frame(stats, row.names = 1)
    sampleinfo <- data.frame(merge(sampleinfo[, setdiff(colnames(sampleinfo), colnames(stats))], stats, by = 0), row.names = 1)

	# update 'history'
	x.hist <- data.frame(history = c("Rarefy:", paste0("Sample size: ", sample.size), 
	                                           paste0("; seed: ",rngseed), 
	                                           paste0("; resampling with replacement: ", replace)), stringsAsFactors = FALSE)
	if (verbose) message("Creating a RepSeqExperiment object...\n")              
	out <- methods::new("RepSeqExperiment", 
	                      assayData = cts, 
	                      sampleData = sampleinfo, 
	                      metaData = mData(x), 
	                      History = rbind(History(x), x.hist))
    cat("Done.\n")
    return(out)
}
