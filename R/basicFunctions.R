utils::globalVariables(c("J", ".", "..sNames", ".SD", ".SDcols", "key", ".N", "count", "..keep.cols", "libnames", "CDR3aa.length", "CDR3aa", "pct", "ctrl.mean", "ID"))

#---------------- basic function ----------------#
#' parse ClonotypeR tables
#'
#' function imports clonotype tables produced by ClonotypeR aligner.
#'
#' @param path full path to a TSV file returned by ClonotypeR (tab-delimited text file, could be gzipped). 
#' @param chain TCR chain \code{A} or [code{B} to extract. Default value is \code{A}.
#' @return a data.table having 7 columns. \code{lib} name of the repertoire, \code{V} V gene identification, \code{J} J gene identification, \code{CDR3aa} CDR3aa chain, \code{CDR3dna} CDR3 DNA chain, \code{score} mapq quality score, \code{count} clonotype assay. Clonotypes were deleted if CDR3aa chain contains STOP codon (*), CDR3dna length is not divisible by 3 or CDR3dna chain contains base "N". 
#' @export
#' @examples
#' \dontrun{
#' dataset <- parseClonotypeR(path, chain="B")
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/rtcr'), 
#'                  package = 'RepSeqData'), 
#'                  full.names = TRUE)
#' l
#' # first element of the vector
#' path <- l[1]
#' dataset <- parseRTCR(path, chain = "B")
#' dataset
#' }
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'              package = 'RepSeqData'), 
#'              full.names = TRUE)
#' l # there are 6 gz-compressed files
#' path <- l[1] # the first gz-compressed of the previous list
#' dataset <- parseMiXCR(path, chain = "B")
#' dataset
#' }
parseMiXCR <- function(path, chain = c("A", "B")) {
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
    data.table::setnames(tab, c("V", "J", "CDR3aa", "CDR3dna", "score", "count"))
    tab[, V := {t1 = gsub("\\*..", "", V); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
    tab[, J := {t1 = gsub("\\*..", "", J); t2 = gsub("\\s*\\([^\\)]+\\)", "", t1); t3 = gsub("\\,.*", "", t2); t4 = gsub("/.*", "", t3)}]
    tab[, CDR3aa := gsub("\\_", "*", CDR3aa)]
    
    ch <- match.arg(chain)
    tab <- switch(ch,
        A = {
            tab[grepl("TRA", V) & grepl("TRA", J), ]
        },
        B = {
            tab[grepl("TRB", V) & grepl("TRB", J), ]    
        })
    sName <- gsub(".tsv|.txt|.gz|.zip|.tar", "", basename(path))
    tab$lib <- sName
    out <- tab[, c("VpJ", "VJ") := list(paste(V, CDR3aa, J), paste(V, J))][, c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score", "count")]
    return(out)
}

#' parse Adaptive tables
#'
#' function imports clonotype tables produced by Adaptive company.
#' 
#' @param path full path to a file returned by Adaptive Biotechnologies.
#' @param chain "A" for alpha or "B" for beta. Default A.
#' @return a data.table 
#' @export
#' @examples
#' \dontrun{
#' adaptive <- parseAdaptive(path)
#' adaptive
#' }
parseAdaptive <- function(path, chain = c("A", "B")) {
    if (path == "" | missing(path))  stop("Empty file name.")
    tab=V=CDR3dna=CDR3aa=V=lib=vMaxResolved=jMaxResolved <- NULL
    if (filetype(path) == "gzfile") {
        tab <- data.table::fread(eval(paste("gunzip -c ", path)))
        } else {
            tab <- data.table::fread(path)
        }
    namelist <- colnames(tab)
    qualCDR3 <- "vIndex"
    aacdr3 <- "aminoAcid"
    dnacdr3 <- "nucleotide"
    keep.cols <- c("vMaxResolved", "jMaxResolved", aacdr3, dnacdr3, qualCDR3, "count (templates/reads)")
    ch <- match.arg(chain)
    tab <- switch(ch,
        A = {
            tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRA", vMaxResolved)]
        },
        B = {
            tab[, ..keep.cols][!(get(aacdr3) == "")][grep("TCRB", vMaxResolved)]
        }
    )
    tab[, vMaxResolved := gsub("\\*..", "", vMaxResolved)][, jMaxResolved := gsub("\\*..", "", jMaxResolved)]
    data.table::setnames(tab, c("V", "J", "CDR3aa", "CDR3dna", "score", "count"))
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
# @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                  package = 'RepSeqData'), 
#'                  full.names = TRUE)
#' l
#' # first element of the vector
#' path <- l[1]
#' dataset <- parseMiXCR(path, chain = "B")
#' dataset
#' datasetFilter <- filterClonotypes(raw = dataset, 
#'                                   keep.ambigurous = FALSE, 
#'                                   keep.ambiguous = FALSE, 
#'                                   keep.unproductive = FALSE, 
#'                                   aa.th = 8, 
#'                                   outFiltered = F)
#' datasetFilter
#' }
filterClonotypes <- function(raw, keep.ambiguous = FALSE, keep.unproductive = FALSE, aa.th = 8, outFiltered = F) {
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
#' @param path a path to a clonotype file, for the moment, only TSV from ClonotypeR is supported. 
#' @param aligner choice between type of input file provided by "Adaptive", "ClonotypeR", "MiXCR", "rTCR". Default is "MiXCR".
#' @param chain a character \code{A} or \code{B} indicating the chain to be imported. Default is \code{A}. 
#' @param keep.ambiguous a boolean choice if ambiguous clonotypes (contain STOP codon) should be kept in analysis. Default is \code{FALSE}. 
#' @param keep.unproductive a boolean choice if unproductive clonotypes (Euclidean dividion of aa length by 3 > 0) should be kept in analysis. Default is \code{FALSE}.
#' @param aa.th an interger indicates the maximal number of amino acids could be deviated from the mean length. Default is \code{8}.
#' @param outFiltered write dropped reads to file. Default folder is getwd().
#' @return a filtered clonotype table.
#' @seealso \code{\link{parseClonotypeR}}, \code{\link{filterClonotypes}}
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                          package = 'RepSeqData'), 
#'                          full.names = TRUE)
#' l
#' # first element of the vector
#' path <- l[1]
#' dataset <- readClonotypes(path = path, 
#'                    aligner = "MiXCR", 
#'                    chain = "B", 
#'                    keep.ambiguous=FALSE, 
#'                    keep.unproductive=FALSE, 
#'                    aa.th=8, 
#'                    outFiltered=F)
#' dataset
#' }
readClonotypes <- function(path, aligner=c("MiXCR", "Adaptive", "rTCR", "ClonotypeR"), chain=c("A", "B"), keep.ambiguous=FALSE, keep.unproductive=FALSE, aa.th=8, outFiltered=F) {
    if (!file.exists(path)) stop("Full path to clonotypes file is required.")
    type <- match.arg(aligner)
    ch <- match.arg(chain)
    raw <- switch(type, 
        ClonotypeR = {
            parseClonotypeR(path, chain = ch) 
        },
        rTCR = { 
            parseRTCR(path)
        },
        MiXCR = {
            parseMiXCR(path, chain = ch)
        },
        Adaptive = {
            parseAdaptive(path, chain = ch)
        })
    out <- filterClonotypes(raw, 
                            keep.ambiguous = keep.ambiguous, 
                            keep.unproductive = keep.unproductive, 
                            aa.th = aa.th, 
                            outFiltered=outFiltered)
    return(out)
}

#
# function assay number of reasd for each clonotype in a repertoire. 
# @param clonotypeTable a data.table of clonotypes returned by function \code{\link{readClonotypes}} which should contain 3 columns naming "V", "J", "pep".
# @return a data.table of assay. 
# @export
# @examples
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
# @examples
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
#' @param aligner software used to align reads ("Adaptive", "rTCR", "ClonotypeR" or "MiXCR"). Default is MiXCR.
#' @param sampleinfo a data frame containing sample information for each clonotype file in fileList. The number of rows of sampleinfo must be identical to the number of file in fileList. A data frame containing If NULL 
#' @param keep.ambiguous a boolean. If TRUE, ambiguous clonotypes (contain STOP codon) will be kept in analysis. Default is \code{FALSE}. 
#' @param keep.unproductive a boolean. If TRUE, unproductive clonotypes (Euclidean dividion of aa length by 3 > 0) will be kept in analysis. Default is \code{FALSE}.
#' @param filterSingleton a boolean. If TRUE clonotypes having 1 count in only 1 sample will be removed. Default is \code{FALSE}.
#' @param aa.th an interger indicates the maximal number of amino acids could be deviated from the mean length. Default is 8.
#' @param raretab a boolean indicating whether a rarefaction tab should be generated. 
#' @seealso \code{\link[RepSeq]{raretabRepSeq}}
#' @return an object of class RepSeqExperiment.
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                      package = 'RepSeqData'), 
#'                      full.names = TRUE)
#' l# list of gz-compressed files 
#' sampleData <- read.table(system.file(file.path('extdata/sampledata.txt'), 
#'                          package='RepSeqData'), 
#'                          sep = "\t", 
#'                          row.names = 2) 
#' dataset <- readClonotypeSet(fileList = l, 
#'                        cores=1L, 
#'                        aligner = "MiXCR", 
#'                        chain = "B", 
#'                        sampleinfo = sampleData, 
#'                        keep.ambiguous=FALSE, 
#'                        keep.unproductive=FALSE, 
#'                        aa.th=8,
#'                        raretab = TRUE)
#' dataset
#' }
readClonotypeSet <- function(fileList, cores = 1L, aligner = c("MiXCR", "Adaptive", "rTCR", "ClonotypeR"), chain = c("A", "B"), 
    sampleinfo = NULL, keep.ambiguous = FALSE, keep.unproductive = FALSE, filterSingleton = FALSE, aa.th = 8, raretab = TRUE) {
    lib=VpJ=V=VJ=V1 <- NULL
    # get number of cores
    cores <- min(parallel::detectCores()-1, cores)    
    cat("Running on", cores, "cores...\n")
    if (length(fileList) == 0) stop("Empty list of files, please check folder path.\n")
    snames <- gsub(".tsv|.txt|.gz|.zip|.tar", "", basename(fileList))
    # selected parser
    parser <- match.arg(aligner)
    # choice of chain to keep
    ch <- match.arg(chain)   
    # parallel run
    if (cores > 1) {   
        Sys.sleep(0.1)
        cat("Loading and filtering clonotypes...\n") 
        cl <- parallel::makeCluster(cores, type = "SOCK", rscript_args = "--vanilla", useXDR = TRUE)
        parallel::clusterExport(cl = cl, varlist = c("filetype", "readClonotypes", "parseClonotypeR", 
                                        "parseRTCR", "parseMiXCR", "parseAdaptive", "filterClonotypes", 
                                        "setnames"))     
        repList <- pbapply::pblapply(cl = cl, 
                                    fileList, 
                                    readClonotypes, 
                                    aligner = parser, 
                                    chain = ch, 
                                    keep.ambiguous = keep.ambiguous, 
                                    keep.unproductive = keep.unproductive, 
                                    aa.th = 8)
        parallel::stopCluster(cl) 
        } else {
            cat("Loading and filtering clonotypes...\n")          
            repList <- lapply(fileList, readClonotypes, aligner = parser, 
                                                        chain = ch, 
                                                        keep.ambiguous = keep.ambiguous, 
                                                        keep.unproductive = keep.unproductive, 
                                                        aa.th = 8)                  
            }
     # count over list 
     cat("Assembling clonotypes...\n")
     countobj <- data.table::rbindlist(repList)
     indx <- which(unlist(lapply(repList, nrow)) == 0)
     if (length(indx) > 0) {
        message(length(indx), " clonotype table(s) have no records after filtering: ", paste0(snames[indx], collapse=", "), ".")
        message("These sample(s) will be excluded.")
     }
     #if (is.null(sampleinfo)) {
	 sdata <- data.frame(Sample = snames, row.names = snames, stringsAsFactors = TRUE)
	 #}
	 if (!is.null(sampleinfo)) sdata <- data.frame(merge(sdata, sampleinfo, by = 0), row.names = 1)
	 if (nrow(sdata) != length(fileList)) stop("Number of files to import differ number of samples in sampleinfo file.")
	 # update columns
	 sdata[sapply(sdata, is.character)] <- lapply(sdata[sapply(sdata, is.character)], as.factor)     
     # get stats about clonotypes
     stats <- data.frame(countobj[, c(.(nReads = sum(count)), lapply(.SD, uniqueN)), .SDcols = c("VpJ", "V", "J", "VJ", "CDR3aa"), by = "lib"], row.names = 1)
     sdata <- data.frame(merge(sdata, stats, by = 0), row.names = 1, stringsAsFactors = TRUE)
     sdata <- sdata[match(rownames(sdata), rownames(stats)),]
     cat("Creating a RepSeqExperiment object...\n") 
	 # Define the 'history'      
	 x.hist <- data.frame(history = c(paste0("data directory=", dirname(fileList[1])), 
	           paste0("readClonotypeSet; cores=", cores, 
	                   "; aligner=", parser, 
	                   "; chain=", ch, 
	                   "; ambiguous ", keep.ambiguous, 
	                   "; unprod ", keep.unproductive,
	                   "; filterSingleton", filterSingleton,
	                   "; aa threshold=", aa.th,
	                   "; raretab", raretab)), 
	                   stringsAsFactors = FALSE)
	 out <- methods::new("RepSeqExperiment", 
	                       assayData = countobj, 
	                       sampleData = sdata, 
	                       metaData = list(), 
	                       History = x.hist)
	 if (raretab == TRUE) {
	   cat("Computing rarefaction table...\n")
	   mData(out)[[1]] <- raretabRepSeq(out)
	   names(mData(out))[1] <- "raretab"
	 }
	 if (filterSingleton) {
	   cat ("Removing singleton clonotypes...")
	   out <- filterassay(out, n = 1)
	 }
     cat("Done.\n")
     return(out)
}

#' build RepSeqExperiment object from clonotype tables
#'
#' function creates RepSeqExperiment object from clonotype tables. 
#' @param clonotypetab a data.table of clonotype tables. The column names should be: lib, V, J, CDR3aa, CDR3dna, VpJ, VJ, score, count.
#' @param sampleinfo a data frame containing sample meta data
#' @return an object of class RepSeqExperiment
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                      package = 'RepSeqData'), 
#'                      full.names = TRUE)
#' l# list of gz-compressed files 
#' sampleData <- read.table(system.file(file.path('extdata/sampledata.txt'), 
#'                      package='RepSeqData'), 
#'                      sep = "\t", 
#'                      row.names = 2) 
#' clontab.tmp <- lapply(l, readClonotypes, aligner="MiXCR")
#' clontab <- data.table::rbindlist(dataset.tmp)
#' dataset <- RepSeqExp(clontab, sampleinfo = sampleData)
#' dataset
#' }
RepSeqExp <- function(clonotypetab, sampleinfo=NULL) {
    coltab <- c("lib", "V", "J", "CDR3aa", "CDR3dna", "VpJ", "VJ", "score", "count")
    if (missing(clonotypetab)) stop("clonotyetable is missing, a clonotype table is expected.")
    if (!is.data.table(clonotypetab)) setDT(clonotypetab)
    if (!all(grepl(paste(coltab, collapse = "|"), colnames(clonotypetab)))) {
        stop("Column names of clonotype table must contain lib, V, J, CDR3aa, CDR3dna, VpJ, VJ, score, count")
        }
    stats <- clonotypetab[, c(.(nReads = sum(count)), lapply(.SD, uniqueN)), by = "lib", .SDcols = c("VpJ", "V", "J", "VJ", "CDR3aa")]
    sNames <- unique(clonotypetab$lib)
    if (is.null(sampleinfo)) {
        sampleinfo <- data.frame(cbind(sample = sNames, stats), row.names = sNames) 
        } else {
            sampleinfo <- data.frame(cbind(sampleinfo, stats), row.names = sNames)
        }
    # setup history    
    x.hist <- data.frame(history = paste0("RepSeqExp; clononotypetab=", 
                        deparse(substitute(clonotypetab)), "; sampleinfo=", 
                        deparse(substitute(sampleinfo))), 
                        stringsAsFactors = FALSE)
    out <- methods::new("RepSeqExperiment", 
                        assayData = clonotypetab, 
                        sampleData = sampleinfo, 
                        metaData = list(), 
                        History = x.hist) 
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), package = 'RepSeqData'), full.names = TRUE)
#' l# list of gz-compressed files 
#' sampleData <- read.table(system.file(file.path('extdata/sampledata.txt'), 
#'                          package='RepSeqData'), 
#'                          sep = "\t", 
#'                          row.names = 2) 
#' dataset <- readClonotypeSet(fileList = l, 
#'                          cores=1L, 
#'                          aligner = "MiXCR", 
#'                          chain = "B",
#'                          sampleinfo = sampleData, 
#'                          keep.ambiguous=FALSE, 
#'                          keep.unproductive=FALSE, 
#'                          aa.th=8)
#' dataset
#' # produce a count matrix of VJ
#' cts <- countFeatures(x = dataset, level = "VJ")
#' cts
#' }
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                    package = 'RepSeqData'), 
#'                    full.names = TRUE)
#' l# list of gz-compressed files 
#' sampleData <- read.table(system.file(file.path('extdata/sampledata.txt'), 
#'                      package='RepSeqData'), 
#'                      sep = "\t", 
#'                      row.names = 2) 
#' dataset <- readClonotypeSet(fileList = l, 
#'                    cores=1L, 
#'                    aligner = "MiXCR", 
#'                    chain = "B", 
#'                    sampleinfo = sampleData, 
#'                    keep.ambiguous=FALSE, 
#'                    keep.unproductive=FALSE, 
#'                    aa.th=8)
#' dataset
#' # produce a count matrix of VJ
#' seg.use <- segmentUsage(x = dataset, level = "VJ")
#' seg.use
#' }
segmentUsage <- function(x, level=c("VpJ", "V", "J", "VJ")) {
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    levelChoice <- match.arg(level)
    cts <- data.table::copy(assay(x))
    cts <- cts[, .(count = sum(count)), by=c("lib", levelChoice)][, prop := count/sum(count), by = "lib"]
    prop <- data.table::dcast(cts, as.formula(paste0(levelChoice, "~lib")), value.var = "prop", fill = 0)
    return(prop)
}

#' filter clonotypes based on counts
#'
#' function filters out every clonotypes having low counts across all samples.
#'
#' @param x an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @param n an integer below this value clonotypes will be filtered out
#' @return an object of class \code{RepSeqExperiment}
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' # filter clonotypes having less than 3 counts accros all samples
#' filterdata <- filterCount(RepSeqData, n=3) 
#' filterdata
#' }
filterCount <- function(x, n=1) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x)) 
    sampleData <- sData(x)
    keep <- cts[, sum(count) <= n, by = VpJ][V1 == FALSE, ]
    res <- cts[keep, on="VpJ"][, V1:=NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), 
    				                    data.frame(history = paste0(nfilter," clonotypes were filtered using filterCount"))), 
    				                    stringsAsFactors=FALSE)
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' # get only singleton clonotypes
#' singletonData <- getSingleton(RepSeqData) 
#' singletonData
#' }
getSingleton <- function(x) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x)) 
    sampleData <- sData(x)
    keep <- cts[, sum(count) == 1, by = VpJ][V1 == TRUE, ]
    if (nrow(keep) == 0) stop("No singletons found.")
    res <- cts[keep, on = "VpJ"][, V1 := NULL]
    setkey(res, lib)
    nfilter <- nrow(unique(cts, by = "VpJ")) - nrow(unique(res, by = "VpJ"))
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), 
    				                    data.frame(history=paste0("function getSingleton: ", nfilter, " singletons"))), 
    				                    stringsAsFactors=FALSE)
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' privateClones <- getPrivates(RepSeqData) # get only private clonotypes
#' provateClones
#' }
getPrivates <- function(x) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    sampleData <- sData(x)
    keep <- cts[, sum(count>0) == 1, by = VpJ][V1 == TRUE]
    res <- cts[keep, on = "VpJ"][, V1 := NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), 
    				                    data.frame(history=paste0(nfilter, " private clonotypes"))), 
    				                    stringsAsFactors=FALSE)
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeq)
#' library(RepSeqData)
#' rownames(sData(RepSeqData))
#' # get overlap clonotypes between the sample S01 and the sample S02 
#' overlapClones <- getOverlaps(RepSeqData, 
#'                            level = "VpJ", 
#'                            libnames=("S01", "S02")) 
#' overlapClones
#' }
getOverlaps <- function(x, level = c("VpJ", "CDR3dna"), libnames = NULL) {
    V1 <- NULL
    if (missing(x)) stop("x is missing.")
    if (!is.RepSeqExperiment(x)) stop("an object of class RepSeqExperiment is expected.")
    cts <- data.table::copy(assay(x))
    sampleData <- sData(x)
    if (is.null(libnames)) libs <- rownames(sampleData) else libs <- libnames 
    keep <- cts[, sum(count > 0) > 1, by = VpJ][V1 == TRUE]
    res <- cts[keep, on = "VpJ"][, V1 := NULL]
    setkey(res, lib)
    nfilter <- nrow(cts) - nrow(res)
    rm(cts, keep)
    sampleData <- sampleData[rownames(sampleData) %in% unique(res$lib), ]
    out <- new("RepSeqExperiment", 
    				assayData = res,
    				sampleData = sampleData,
    				metaData = mData(x),
    				History = data.frame(rbind(History(x), 
    				            data.frame(history = paste0(nfilter, " shared clonotypes"))), 
    				            stringsAsFactors=FALSE)
	   	   )
    return(out)    
}

#' filter clonotypes presented in less than a proportion of samples.
#'
#' A function that returns a function with values for A, p and na.rm bound to the specified values. The function takes a single vector, x, as an argument. When the returned function is evaluated it returns TRUE if the proportion of values in x that are larger than A is at least p.
#'
#' @param x an object of class RepSeqExperiment.
#' @param freq a number between 0 & 1. A clonotype is retained if it is present in a proportion of samples that excces This is the proportion of samples 
#' @param group a factor in sampleData, if null, freq is applied to the group having the smallest number of membres.
#' @return an object of class RepSeqExperiment satisfying the conditions.
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeq)
#' library(RepSeqData)
#' sData(RepSeqData)
#' # get overlap clonotypes between the sample S01 and the sample S02
#' filterdata <- filterFrequency(RepSeqData, 
#'                        freq = 0.2, 
#'                        group = "quantity")  
#' filterdata
#' }
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
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' l <- list.files(system.file(file.path('extdata/mixcr'), 
#'                          package = 'RepSeqData'), 
#'                          full.names = TRUE)
#' l
#' # first element of the vector
#' filelist1 <- l[1:3]
#' dataset1 <- readClonotypeSet(fileList = filelist1, 
#'                    aligner = "MiXCR", 
#'                    chain = "B", 
#'                    keep.ambiguous=FALSE, 
#'                    keep.unproductive=FALSE, 
#'                    aa.th=8)
#' dataset1
#' filelist2 <- l[4:6]
#' dataset2 <- readClonotypeSet(fileList = filelist2, 
#'                    aligner = "MiXCR", 
#'                    chain = "B", 
#'                    keep.ambiguous=FALSE, 
#'                    keep.unproductive=FALSE, 
#'                    aa.th=8)
#' dataset2
#' dataset <- concateRepSeq(a = dataset1, b = dataset2)
#' dataset
#' }
concateRepSeq <- function(a, b) {
    if (missing(a) | missing (b)) stop("Two RepSeqExperiment objects are required.")
    if (!is.RepSeqExperiment(a)) stop("a is not an object of class RepSeqExperiment.")
    if (!is.RepSeqExperiment(b)) stop("b is not an object of class RepSeqExperiment.")
    if (any(rownames(RepSeq::sData(a)) == rownames(RepSeq::sData(b)))) stop("Common sample names are not allowed, please use the function names()<- to rename them.")
    cts <- rbind(RepSeq::assay(a), RepSeq::assay(b))
    cts[, lib := as.character(lib)]
    sampleinfo <- rbind(RepSeq::sData(a), RepSeq::sData(b))
    a.history <- paste0(deparse(substitute(a)), ":", RepSeq::History(a)$history)
    b.history <- paste0(deparse(substitute(b)), ":", RepSeq::History(b)$history)
    concat.history <- paste(date(),"- concatenation of", deparse(substitute(a)), "and", deparse(substitute(b)), "using the function concateRepSeq")
    all.history <- data.frame(history=c(a.history, b.history, concat.history))
    metainfo <- ifelse(length(RepSeq::mData(a)) > 0 | length(RepSeq::mData(b)) > 0, c(RepSeq::mData(a), RepSeq::mData(b)), list())
    out <- new("RepSeqExperiment", 
            assayData = cts, 
            sampleData = sampleinfo, 
            metaData = metainfo, 
            History = all.history)
    return(out)
}

#' drop sample(s) from a RepSeqExperiment object
#'
#' function allows to remove one or several samples from an RepSeqExperiment object. 
#' @param x a RepSeqExperiment object.
#' @param samples a vector containing sample names to be removed or a vector integer of their positions in sampleData slot.
#' @return a RepSeqExperiment object.
#' @export
#' @examples
#' \dontrun{
#' # The package RepSeqData contains example datasets 
#' library(RepSeqData)
#' # show sample names
#' rownames(sData(RepSeqData)
#' subsetRepSeqData <- dropSamples(x = RepSeqData, samples=c("S10", "S11", "S12")) 
#' subsetRepSeqData
#' }
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
    cts[, lib := as.character(lib)]
    sampleinfo <- droplevels(sampleinfo[-c(index), , drop=FALSE])
    x.history <- data.frame(rbind(History(x), data.frame(history = paste(date(), "- drop", paste0(snames, collapse=", "), "from", 
                                                            deparse(substitute(x)), "using the function dropSamples"))))
    metainfo <- mData(x)
    out <- new("RepSeqExperiment", 
            assayData = cts, 
            sampleData = sampleinfo, 
            metaData = metainfo, 
            History = x.history)
    return(out)
}

