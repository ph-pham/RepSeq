#-------- diversity computation -----------#

#' updates sample data with basic statistics
#'
#' function produces several useful diversity indices. 
#'
#' @param x an object of class RepSeqExperiment
#' @param level level of reppertoire to assess diversity, choice are VpJ, V, J, VJ or CDR3aa)
#' @return a data frame of Chao
#' @export
basicIndices <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    # summary function
    pastek <- function(y) list(shannon=.diversity(y, index="shannon"), simpson=.diversity(y, index="simpson"), 
                            invsimpson=.diversity(y, index="invsimpson"), gini=.gini(y), 
                            chao1=.chao1(y)$chao.est, chao1.se=.chao1(y)$chao.se, ichao=iChao(y), chaowor=Chaowor(y))   
    out <- assay(x)[, .(count=sum(count)), by=c("lib", levelChoice)][, pastek(count), by="lib"]
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
divLevel <- function(x, index=c("shannon", "simpson", "invsimpson"), level=c("VJ", "V", "J"), norm=TRUE, base=exp(1)) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
    divIndex <- match.arg(index)
    levelChoice <- match.arg(level)
    res <- assay(x)[, .diversity(count, index=index, base=base, norm=norm), by=c("lib", levelChoice)]
    out <- data.table::dcast(res, as.formula(paste0(levelChoice, "~lib")), value.var="V1", fill=0)
    return(out)
}

# compute renyi or hill number
#
# function computes Renyi (Hill's numbers) indices according to \alpha.
# @param x a vector of count
# 
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
# @example
# renyiProfiles()
renyiProfiles <- function(x, scales=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), level=c("VpJ", "V", "J", "VJ", "CDR3aa")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is exprected.")
    levelChoice <- match.arg(level)
    out <- assay(x)[, .(count=sum(count)), by=c("lib", levelChoice)][, lapply(scales, function(y) .renyiCal(count, y)), by="lib"]
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
# @param x an object of class RepSeqExperiment
#
# 
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

#.d50 <- function(x) {
#    
#
#}

#' Improved Chao1
#'
#' function computes the improve version of Chao1
#' @param x a vector of count.
#' @return improved Chao1 value (ref) Chao and Lin Chao, A. and Lin, C.-W. (2012) Nonparametric lower bounds for species richness and shared species richness under sampling without replacement. Biometrics,68, 912–921. 
#' @export
# @example
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
# @example
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


