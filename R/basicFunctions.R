utils::globalVariables(c("J", ".", "..sNames", ".SD", ".SDcols", "key", ".N", "count", "..keep.cols", "libnames"))

#---------------- basic function ----------------#
#' parse ClonotypeR tables
#'
#' function imports clonotype tables produced by ClonotypeR aligner.
#'
#' @param path full path to a TSV file returned by ClonotypeR (tab-delimited text file, could be gzipped). 
#' @param chain TCR chain \code{A} or [code{B} to extract. Default value is \code{A}.
#' @return a data.table having 7 columns. \code{lib} name of the repertoire, \code{V} V gene identification, \code{J} J gene identification, \code{CDR3aa} CDR3aa chain, \code{CDR3dna} CDR3 DNA chain, \code{score} mapq quality score, \code{count} clonotype assay. Clonotypes were deleted if CDR3aa chain contains STOP codon (*), CDR3dna length is not divisible by 3 or CDR3dna chain contains base "N". 
#' @export
#' @example
#' \dontrun{
#' parseClonotypeR()
#' }
parseClonotypeR <- function(path, chain=c("A", "B")) {
    if (path == "" | missing(path))  stop("Empty file name.")
    tab=V=CDR3dna=CDR3aa=V=lib <- NULL
    if (filetype(path) == "gzfile") {
        tab <- data.table::fread(eval(paste("gunzip -c ", path)))
        } else {
            tab <- data.table::fread(path)
        }
    data.table::setnames(tab, c("lib", "V", "J", "score", "mapq", "read", "CDR3dna", "Quality", "CDR3aa"))
    ch <- match.arg(chain)
    if (ch == "A") {
        tab <- tab[grepl("TRA", V) & grepl("TRA", J), ]
        }
    if (ch == "B") {
        tab <- tab[grepl("TRB", V) & grepl("TRB", J), ]    
        }
    tab[, c("lib", "VpJ", "VJ") := list(gsub(".tsv|.txt|.gz|.zip|.tar", "", lib), paste(V, CDR3aa, J), paste(V, J))]   
    out <- tab[, c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score") ]
    #data.table::setkey(out, lib, V, J, CDR3aa, CDR3dna, VpJ, VJ, score)
    out <- out[, .(count = .N), by = c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score")]
    #data.table::setnames(out, c("lib", "V", "J", "CDR3aa", "VpJ", "VJ", "score", "count"))
    return(out)
}


#' parse rTCR tables
#'
#' function imports clonotype tables produced by rTCR aligner.
#'
#' @param path full path to a file returned by rTCR
#' @return a data.table 
#' @export
# @example
# run parseRTCR()
parseRTCR <- function(path) {
    if (missing(path)) stop("path to rtcr output is missing.")
    V=CDR3aa <- NULL
    # load     
    if (filetype(path) == "gzfile") {
        rtcrtab <- data.table::fread(eval(paste("gunzip -c ", path)))
        } else {
            rtcrtab <- data.table::fread(path)
        }
    sName <- gsub(".tsv|.txt|.gz|.zip|.tar", "", basename(path))
    data.table::setnames(rtcrtab, c("count", "CDR3aa", "V", "J", "CDR3dna", "V.pos", "J.pos", "Frame", "stopCodons", "score", "Quality"))
    rtcrtab$lib <- sName
    out <- rtcrtab[, c("V", "J") := list(gsub("\\*..", "", V), gsub("\\*..", "", J))][, c("VpJ", "VJ") := list(paste(V, CDR3aa, J), paste(V, J))]
    out <- out[, c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score", "count")]
    return(out)
}

#' parse MiXCR tables
#'
#' function imports clonotype tables produced by MiXCR aligner.
#' 
#' @param path full path to a file returned by MiXCR
#' @param chain TCR chain \code{A} or [code{B} to extract. Default value is \code{A}.
#' @return a data.table 
#' @export
# @example
# run parseMiXCR()
parseMiXCR <- function(path, chain=c("A", "B")) {
    if (path == "" | missing(path))  stop("Empty file name.") 
    tab=V=CDR3dna=CDR3aa=V=lib <- NULL
    if (filetype(path) == "gzfile") {
        tab <- data.table::fread(eval(paste("gunzip -c ", path)))
        } else {
            tab <- data.table::fread(path)
        }
    namelist <- colnames(tab)
    qualCDR3 <- ifelse(length(grep("minQualShortCDR3", namelist)) > 0, "minQualShortCDR3", "minQualCDR3")
    aacdr3 <- ifelse(length(grep("aaSeqShortCDR3", namelist)) > 0, "aaSeqShortCDR3", "aaSeqCDR3")
    dnacdr3 <- ifelse(length(grep("nSeqShortCDR3", namelist)) > 0, "nSeqShortCDR3", "nSeqCDR3")
    keep.cols <- c("allVHitsWithScore", "allJHitsWithScore", aacdr3, dnacdr3, qualCDR3, "cloneCount")
  
    tab <- tab[, ..keep.cols]
    #tab <- tab[, c("allVHitsWithScore", "allJHitsWithScore", "aaSeqCDR3", "nSeqCDR3", "minQualCDR3", "cloneCount")]
    data.table::setnames(tab, c("V", "J", "CDR3aa", "CDR3dna", "score", "count"))
    tab[, V := {t1 = gsub("\\*..", "", V); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
    tab[, J := {t1 = gsub("\\*..", "", J); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
    tab[, CDR3aa := gsub("\\_", "*", CDR3aa)]
    
    ch <- match.arg(chain)
    if (ch == "A") {
        tab <- tab[grepl("TRA", V) & grepl("TRA", J), ]
        }
    if (ch == "B") {
        tab <- tab[grepl("TRB", V) & grepl("TRB", J), ]    
        }
    sName <- gsub(".tsv|.txt|.gz|.zip|.tar", "", basename(path))
    tab$lib <- sName
    out <- tab[, c("VpJ", "VJ") := list(paste(V, CDR3aa, J), paste(V, J))][, c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score", "count")]
    return(out)
}


#' filter clonotypes
#'
#' function loads a clonotype table and applies a filter according to several conditions. 
#'
#' @param raw a clonotype table.
#' @param keep.ambiguous a boolean choice if ambiguous clonotypes (contain STOP codon) should be kept in analysis. Default is \code{FALSE}. 
#' @param keep.unproductive a boolean choice if unproductive clonotypes (Euclidean dividion of aa length by 3 > 0) should be kept in analysis. Default is \code{FALSE}.
#' @param aa.th an interger indicates the maximal number of amino acids could be deviated from the mean length. Default is \code{8}.
#' @param outFiltered write dropped reads to file. Default folder is getwd().
#' @return a filtered clonotype table.
#' @export
# @example
# run filterClonotypes()
filterClonotypes <- function(raw, keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=8, outFiltered=F) {
    CDR3aa=CDR3dna <- NULL
    #cat(paste("Processing", x, "\n"))
    if (missing(raw)) stop("a clonotype table is required.\n")
    if (!data.table::is.data.table(raw)) stop("a data.table is expected.\n")
    # distribution of aa length
    aa.size <- raw[, nchar(CDR3aa)]
    # frequency of aa lengths
    aa.distr <- table(aa.size)
    # get aa size 
    aa.length <- as.numeric(names(aa.distr))
    aa.mean <- aa.length[which.max(aa.distr)]
    # keep read
    keep.length <- aa.size <= (aa.mean + aa.th)
    # read to keep
    indx <- keep.length
    if (!keep.ambiguous)
        indx <- indx & !raw[, grepl("N", CDR3dna)]
    if (!keep.unproductive)
        indx <- indx & !raw[, nchar(CDR3dna) %% 3 > 0 | grepl("\\*", CDR3aa)]
    # write dropped reads to file 
    if (outFiltered) data.table::fwrite(raw[!indx, ], file=file.path(getwd(), paste("drop_", raw$lib[1], ".csv", sep="")), row.names=TRUE) 
    # return
    out <- raw[indx, ]
    return(out)
}

#' read clonotype table
#'
#' function readColontypes is a wrapper of parseClonotypes & filterClonotypes
#'
#' @param x a path to a clonotype file, for the moment, only TSV from ClonotypeR is supported. 
#' @param aligner choice between type of input file provided by "ClonotypeR", "MiXCR", "rTCR". Default is "ClonotypeR".
#' @param chain a character \code{A} or \code{B} indicating the chain to be imported. Default is \code{A}. 
#' @param keep.ambiguous a boolean choice if ambiguous clonotypes (contain STOP codon) should be kept in analysis. Default is \code{FALSE}. 
#' @param keep.unproductive a boolean choice if unproductive clonotypes (Euclidean dividion of aa length by 3 > 0) should be kept in analysis. Default is \code{FALSE}.
#' @param aa.th an interger indicates the maximal number of amino acids could be deviated from the mean length. Default is \code{8}.
#' @param outFiltered write dropped reads to file. Default folder is getwd().
#' @return a filtered clonotype table.
#' @seealso \code{\link{parseClonotypeR}}, \code{\link{filterClonotypes}}
#' @export
# @example
# run readClonotypes()
readClonotypes <- function(x, aligner=c("ClonotypeR", "rTCR", "MiXCR"), chain=c("A", "B"), keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=8, outFiltered=F) {
    if (!file.exists(x)) stop("Full path to clonotypes file is required.")
    type <- match.arg(aligner)
    ch <- match.arg(chain)
    if (type == "ClonotypeR") raw <- parseClonotypeR(x, chain=ch)
    if (type == "rTCR") raw <- parseRTCR(x)
    if (type == "MiXCR") raw <- parseMiXCR(x, chain=ch)   
    out <- filterClonotypes(raw, keep.ambiguous=keep.ambiguous, keep.unproductive=keep.unproductive, aa.th=aa.th, outFiltered=outFiltered)
    return(out)
}

#
# function assay number of reasd for each clonotype in a repertoire. 
# @param clonotypeTable a data.table of clonotypes returned by function \code{\link{readClonotypes}} which should contain 3 columns naming "V", "J", "pep".
# @return a data.table of assay. 
# @export
# @example
# readClonotypes()
# run countVpJ()
#countVpJ <- function(clonotypeTable) {
#    if (missing(clonotypeTable)) stop("Clonotype table is required.\n")
#    if (!data.table::is.data.table(clonotypeTable)) stop("a data.table is required.")
#    keys <- c("V", "pep", "J")
#    if (sum(keys %in% colnames(clonotypeTable)) != length(keys)) stop("ClonotypeTable musts contain columns named V, pep and J.")
#    col2drop <- setdiff(colnames(clonotypeTable), keys)
#    V=VpJ=pep=J <- NULL 
#    if (length(col2drop>0)) {
#        clonotypeTable[, (col2drop) := NULL]
#        }
#    clonotypeTable[, VpJ := paste(V, pep, J)]
#    res <- clonotypeTable[, .(count=.N), by = VpJ]
#    data.table::setkey(res, VpJ)
#    return(res)
#}

# function mergeVpJcount assay clonotypes for several repertoires and merges them together.
# @param repList a list of clonotype tables
# @return a marix of assay with \code{V-pep-J} in rows and samples in columns.
# @seealso \code{\link{countVpJ}}, \code{\link{pbapply}}
# @export
# @example
# run mergeVpJcount()
#mergeVpJcount <- function(countList) {
#    if (missing(countList)) stop("A list of clonotype tables is required.\n")
#    # initialization
#    VpJ=V=pep=J=VJ <- NULL
#    # get unique list of V-pep-J
#    lib <- unique(do.call('c', lapply(countList, function(x) x[[1]] )))
#    # list of assay according to the lib (union of V-pep-J), return count if matched otherwise NA
#    count.list <- lapply(countList, function(x) `[[`(x,2)[data.table::chmatch(lib, `[[`(x,1))])
#    # count Matrix
#    # Merging assay
#    countMatrix <- do.call('cbind', count.list)
#    # set names
#    sNames <- names(countList)
#    rownames(countMatrix) <- lib
#    colnames(countMatrix) <- sNames
#    # na to 0
#    countMatrix[is.na(countMatrix)] <- 0
#    countMatrix <- data.table::data.table(countMatrix, keep.rownames="VpJ", key="VpJ")
#    countMatrix[, c("V", "pep", "J"):=tstrsplit(VpJ, " ", fixed=TRUE)][, VJ:=paste(V, J)]
#    data.table::setkey(countMatrix, VpJ, VJ, V, pep, J)   
#    data.table::setcolorder(countMatrix, c(data.table::key(countMatrix),sNames))
#    return(countMatrix)
#}

#' read clonotype table set
#'
#' function reads all clonotype tables stored in a folder, filters reads and assay reads for each clonotype.
#'
#' @param fileList a list of file paths to import. 
#' @param cores an interger number of cores to use. Default is 1.
#' @param chain an character indicates which TCR chain alpha or beta to import, Use \code{A} for alpha chain and \code{B} for beta chain. Default is \code{A}.
#' @param aligner software used to align reads ("rTCR", "ClonotypeR" or "MiXCR"). 
#' @param sampleinfo a data frame containing sample information for each clonotype file in fileList. The number of rows of sampleinfo must be identical to the number of file in fileList. A data frame containing If NULL 
#' @param keep.ambiguous a boolean choice if ambiguous clonotypes (contain STOP codon) should be kept in analysis. Default is \code{FALSE}. 
#' @param keep.unproductive a boolean choice if unproductive clonotypes (Euclidean dividion of aa length by 3 > 0) should be kept in analysis. Default is \code{FALSE}.
#' @param aa.th an interger indicates the maximal number of amino acids could be deviated from the mean length. Default is 8.
#' @return an object of class RepSeqExperiment.
#' @export
# @example
# readClonotypeSet()
readClonotypeSet <- function(fileList, cores=1L, aligner=c("ClonotypeR", "rTCR", "MiXCR"), chain=c("A", "B"), 
    sampleinfo=NULL, keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=8) {
    lib=VpJ=V=VJ <- NULL
    # get number of cores
    cores <- min(parallel::detectCores()-1, cores)    
    cat("Running on", cores, "cores.\n")
    if (length(fileList) == 0) stop("Empty list of files, please check folder path.\n") 
    if (is.null(sampleinfo)) {
	   sampleinfo <- data.frame(ID=basename(fileList), row.names=gsub(".tsv|.txt|.gz|.zip|.tar", "", basename(fileList)))
	   }
	if (nrow(sampleinfo) != length(fileList)) stop("Number of files to import differ number of samples in sampleinfo file.")
    # selected parser
    parser <- match.arg(aligner)
    # choice of chain to keep
    ch <- match.arg(chain)   
    # parallel run
    if (cores > 1) {   
        Sys.sleep(0.1)
        cat("Loading and filtering clonotypes...\n") 
        cl <- parallel::makeCluster(cores, setup_strategy = "sequential")
        parallel::clusterExport(cl=cl, varlist=c("filetype", "readClonotypes", "parseClonotypeR", "parseRTCR", "parseMiXCR", "filterClonotypes", "setnames"))     
        repList <- pbapply::pblapply(cl=cl, fileList, readClonotypes, aligner=parser, chain=ch, keep.ambiguous=keep.ambiguous, keep.unproductive=keep.unproductive, aa.th=8)
        parallel::stopCluster(cl) 
        } else {
            cat("Loading and filtering clonotypes...\n")          
            repList <- lapply(fileList, readClonotypes, aligner=parser, chain=ch, keep.ambiguous=keep.ambiguous, keep.unproductive=keep.unproductive, aa.th=8)
            #names(repList) <- basename(fileList)                  
            }
     # count over list 
     cat("Assembling clonotypes...\n")
     countobj <- data.table::rbindlist(repList)
     # get stats about clonotypes
     stats <- countobj[, c(.(nReads=sum(count)), lapply(.SD, uniqueN)), by="lib", .SDcols=c("VpJ", "V", "J", "VJ", "CDR3aa")]
     sampleinfo <- data.frame(sampleinfo, data.frame(stats, row.names=1))
     # make keys
     #data.table::setkey(countobj, lib, VpJ, V, J, VJ)
     cat("Creating a RepSeqExperiment object...\n") 
	 	 # Define the 'history'      
	 x.hist <- data.frame(history=c(paste0("data directory=", dirname(fileList[1])), 
	           paste0("readClonotypeSet; cores=", cores, "; aligner=", parser, "; chain=", ch, "; ambiguous ", 
	               keep.ambiguous, "; unprod ", keep.unproductive, "; aa threshold=", aa.th)), stringsAsFactors=FALSE)
	 out <- methods::new("RepSeqExperiment", assayData=countobj, sampleData=sampleinfo, metaData=list(), History=x.hist)
     cat("Done.\n")
     return(out)
}

# create V, pep, J columns
#
# function makeFeatures creates feature data from row names of count data (under the format V-p-J) the column VpJ of count data.
#
# @param x a character vector, a matrix, data frame or data.table having row names under the format "V pep J".
# @return a data.table containing features.
# @export
# @example
# makeFeatures()
#makeFeatures <- function(x) {
#    if (missing(x)) stop ("x is missing.")
#    # initialization
#    VpJ = V = pep = J = VJ <- NULL
#    if (is.matrix(x) || is.data.frame(x)) {
#        features <- data.table::data.table(VpJ=rownames(x))
#   }
#   if (data.table::is.data.table(x)) {
#        if (is.null(x$VpJ)) stop("data.table x must contain a column named ")
#        features <- data.table::data.table(VpJ=x$VpJ)
#   }
#   features <- features[, c("V", "pep", "J") := data.table::tstrsplit(VpJ, " ", fixed=TRUE)]
#   features <- features[, VJ := paste(V, J)]
#   return(features)
#}


#' build RepSeqExperiment object from clonotype tables
#'
#' function creates RepSeqExperiment object from clonotype tables. 
#' @param clonotypetab a data.table of clonotype tables. The column names should be: lib, V, J, CDR3aa, CDR3dna, VpJ, VJ, score, count.
#' @param sampleinfo a data frame containing sample meta data
#' @return an object of class RepSeqExperiment
#' @export
# @example
RepSeqExp <- function(clonotypetab, sampleinfo=NULL) {
    coltab <- c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score", "count")
    if (missing(clonotypetab)) stop("clonotyetable is missing, a clonotype table is expected.")
    if (!is.data.table(clonotypetab)) setDT(clonotypetab)
    if (!all(grepl(paste(coltab, collapse="|"), colnames(clonotypetab)))) {
        stop("Column names of clonotype table must contain lib, V, J, CDR3aa, CDR3dna, VpJ, VJ, score, count")
        }
    stats <- clonotypetab[, c(.(nReads=sum(count)), lapply(.SD, uniqueN)), by="lib", .SDcols=c("VpJ", "V", "J", "VJ", "CDR3aa")]
    sNames <- unique(clonotypetab$lib)
    if (is.null(sampleinfo)) {
        sampleinfo <- data.frame(cbind(sample=sNames, stats), row.names=sNames) 
        } else {
            sampleinfo <- data.frame(cbind(sampleinfo, stats), row.names=sNames)
        }
    x.hist <- data.frame(history = paste0("RepSeqExp; clononotypetab=", deparse(substitute(clonotypetab)), "; sampleinfo=", deparse(substitute(sampleinfo))), stringsAsFactors=FALSE)
    out <- methods::new("RepSeqExperiment", assayData=clonotypetab, sampleData=sampleinfo, metaData=list(), History=x.hist) 
    return(out)
}

#' count reads
#'
#' function assay reads by V, J or V-J genes.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param level "V", "J" or "VJ" genes.
#' @return a data.table of assay with unique features in rows and samples in columns.
#' @export
# @example
# run countFeatures()
countFeatures <- function(x, level=c("VpJ", "V", "J", "VJ", "CDR3aa")) {    
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("An object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    out <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~ lib")), value.var="count", fun=sum)
    return(out)
}

#' compute usage of segments
#'
#' This function computed clonotype usage 
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}], returned by \code{\link{readClonotypeSet}}
#' @param level segment usage (ie percentage) of VpJ, V, J, or VJ within a repertoire.
#' @return a data.table of segment usage of genes in each repertoire.
#' @export
# @example
# segmentUsage()
segmentUsage <- function(x, level=c("VpJ", "V", "J", "VJ")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    cts <- cts[, .(count=sum(count)), by=c("lib", levelChoice)][, prop:=count/sum(count), by="lib"]
    prop <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~lib")), value.var="prop", fill=0)
    return(prop)
}

#' filter clonotypes based on assay
#'
#' function filters out every clonotypes having low assay across all samples.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param n an integer below this value clonotypes will be filtered out
#' @return an object of class \code{RepSeqExperiment}
#' @export
# @example
# filterassay
filterassay <- function(x, n=1) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x)) 
    sampleData <- sData(x)
    keep <- cts[, sum(count) <= n, by=VpJ][V1 == FALSE, ]
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), data.frame(history=paste0(nfilter," clonotypes were filtered using filterassay"))), stringsAsFactors=FALSE)
	   			)
    out
}

#' get singletons 
#'
#' function gets all clonotype having only 1 count across all samples
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @return an object of class \code{RepSeqExperiment}
#' @export
# @example
# getSingleton
getSingleton <- function(x) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x)) 
    sampleData <- sData(x)
    keep <- cts[, sum(count) == 1, by=VpJ][V1 == TRUE, ]
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), data.frame(history=paste0(nfilter, " singletons"))), stringsAsFactors=FALSE)
	   			)
    return(out)
}


#' get private clonotypes
#' 
#' function returns private clonotyes which expressed only in one sample
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @return an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @export
# @example
getPrivates <- function(x) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    sampleData <- sData(x)
    keep <- cts[, sum(count>0)==1, by=VpJ][V1==TRUE, ]
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), data.frame(history=paste0(nfilter, " private clonotypes"))), stringsAsFactors=FALSE)
	   			)
    return(out)    
}

#' get shared clonotypes
#' 
#' function returns an RepSeqExperiment object containing shared clonotyes which expressed in at least two samples
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param level level of shared clonotypes, VpJ or CDR3dna  
#' @param libnames a vector of specific sample names to get shared clonotypes, default value is NULL, shared clonotypes will be computed for all samples.
#' @return an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @export
# @example
getOverlaps <- function(x, level=c("VpJ", "CDR3dna"), libnames=NULL) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    sampleData <- sData(x)
    if (is.null(libnames)) libs <- rownames(sampleData) else libs <- libnames 
    keep <- cts[, sum(count>0)>1, by=VpJ][V1==TRUE, ]
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), data.frame(history=paste0(nfilter, " shared clonotypes"))), stringsAsFactors=FALSE)
	   			     )
    return(out)    
}


#' filter based on expression.
#'
#' this function allows to filter data according to the percentage of presence across samples.
#'
#' @param x count data set with features in rows and samples in columns.
#' @param freq a threshold expressed between 0 & 1 below which data will be filtered out.
#' @param group a factor in sampleData, if not null, freq is applied to the group having the smallest number of membres. 
#' @return subset of x that match the conditions.
#' @export
# @example
# filterFrequency()
filterFrequency <- function(x, freq=0.2, group=NULL) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    if (freq <= 0 | freq > 1) stop("freq musts be a number between 0 and 1.")   
    sampleData <- sData(x)
    cts <- data.table::copy(assay(x))
    if (!is.null(group)) { 
        grp <- sampleData[, group]
        n <- min(table(grp))
        grp.min.name <- names(which.min(table(grp)))
        libnames <- rownames(sampleData[grp %in% grp.min.name, ])
        keep <- cts[lib %in% libnames, sum(count>0) > n * freq, by=VpJ][V1==TRUE, ]
        x.hist <- paste0("filterFrequency; freq=", freq, "; group=", group, "; (smallest group: ", grp.min.name, ")")
        } else {
            n <- nrow(sampleData)
            keep <- cts[, sum(count>0) > n * freq, by=VpJ][V1==TRUE, ]
            x.hist <- paste0("filterFrequency; freq=", freq, "; group=", group)
            }
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment",
        assayData = res,
        sampleData = sampleData,
        metaData = mData(x),
        History = rbind(History(x), data.frame(history=x.hist))
        )
    return(out)
}

#' get type of file
#'
#' function allows to get the type of a file
#' @param path path to a file
#' @return a string indicating file type
#' @export
filetype <- function(path) {
    f <- file(path)
    ext <- summary(f)$class
    close.connection(f)
    return(ext)
}

#' concatenate 2 RepSeqExperiment objects
#'
#' function concateRepSeq allow to concatenate 2 RepSeqExperiement objects. Overlaps in sample names are not allowed. 
#' @param a first RepSeqExperiment object.
#' @param b second RepSeqExperiment object.
#' @return a RepSeqExperiment object.
#' @export
# @example
concateRepSeq <- function(a, b) {
    if (missing(a) | missing (b)) stop("Two RepSeqExperiment objects are required.")
    if (!is.RepSeqExperiment(a)) stop("a is not an object of class RepSeqExperiment.")
    if (!is.RepSeqExperiment(b)) stop("b is not an object of class RepSeqExperiment.")
    if (any(rownames(RepSeq::sData(a)) == rownames(RepSeq::sData(b)))) stop("Common sample names are not allowed, please use the function names()<- to rename them.")
    cts <- rbind(RepSeq::assay(a), RepSeq::assay(b))
    cts[, lib:=as.character(lib)]
    sampleinfo <- rbind(RepSeq::sData(a), RepSeq::sData(b))
    a.history <- paste0(deparse(substitute(a)), ":", RepSeq::History(a)$history)
    b.history <- paste0(deparse(substitute(b)), ":", RepSeq::History(b)$history)
    concat.history <- paste(date(),"- concatenation of", deparse(substitute(a)), "and", deparse(substitute(b)), "using the function concateRepSeq")
    all.history <- data.frame(history=c(a.history, b.history, concat.history))
    metainfo <- ifelse(length(RepSeq::mData(a)) > 0 | length(RepSeq::mData(b)) > 0, c(RepSeq::mData(a), RepSeq::mData(b)), list())
    out <- new("RepSeqExperiment", 
            assayData=cts, 
            sampleData=sampleinfo, 
            metaData=metainfo, 
            History=all.history)
    return(out)
}

#' drop sample(s) from a RepSeqExperiment object
#'
#' function allows to remove one or several samples from an RepSeqExperiment object. 
#' @param x a RepSeqExperiment object.
#' @param samples a vector containing sample names to be removed or a vector of their positions in sampleData slot.
#' @return a RepSeqExperiment object.
#' @export
# @example
dropSamples <- function(x, samples) {
    if (missing(x)) stop("A RepSeqExperiment object is required.")
    if (!is.RepSeqExperiment(x)) stop("x is not an object of class RepSeqExperiment.")
    sampleinfo <- sData(x)
    if (is.numeric(samples)) {
        index <- samples
        snames <- rownames(sampleinfo)[index]
    }
    if (is.character(samples)) {
        if (!all(samples %in% rownames(sampleinfo))) stop("Sample names not found in x.")
        index <- which(rownames(sampleinfo) %in% samples)
        snames <- samples
    }
    cts <- copy(assay(x))
    cts <- cts[!(lib %in% snames)]
    cts[, lib:=as.character(lib)]
    sampleinfo <- droplevels(sampleinfo[-c(index), , drop=FALSE])
    x.history <- data.frame(rbind(History(x), data.frame(history=paste(date(), "- drop", paste0(snames, collapse=", "), "from", deparse(substitute(x)), "using the function dropSamples"))))
    metainfo <- mData(x)
    out <- new("RepSeqExperiment", 
            assayData=cts, 
            sampleData=sampleinfo, 
            metaData=metainfo, 
            History=x.history)
    return(out)
}

