# without this, there is a warning during check
# data.table variables
utils::globalVariables(c("J", ".", "VpJ", "lib", "V", "J", "VJ", "..cols", "..n.."))

#------------------------------------------------------------------
# Define the class
#------------------------------------------------------------------
#' Class RepSeqExperiment
#'
#' An S4 class to represent a HTS RepSeq set.
#'
#' @slot assayData a data.table containing information about clonotype.
#' @slot sampleData a data frame containing sample information such as treatment groups, species, ...
#' @slot metaData a list of meta data
#' @slot History a data frame registering operations done on the object.
#' @author H. P. Pham (mailto: hp.pham(a)iltoopharma.fr, twitter: #kiniou.est.mon.ami, #chokomlabreizh)
#' @rdname RepSeqExperiment-class
#' @name RepSeqExperiment-class
#' @exportClass RepSeqExperiment
RepSeqExperiment <- setClass("RepSeqExperiment", 
    representation = representation(
        assayData = "data.table", 
        sampleData = "data.frame",
        metaData = "list",
        History = "data.frame"
        )
)

#' Method assay.
#'
#' @param object a RepSeqExperiment object
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod assay
setGeneric("assay", function(object) standardGeneric("assay"))

#' Method assay<-.
#'
#' @param value either numeric, character or data frame 
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod assay<-
setGeneric("assay<-", function(object, i, j, value) standardGeneric("assay<-"))

#' Method mData.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod mData
setGeneric("mData", function(object) standardGeneric("mData"))

#' Method mData<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod mData<-
setGeneric("mData<-", function(object, value) standardGeneric("mData<-"))

#' Method sData.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod sData
setGeneric("sData", function(object) standardGeneric("sData"))

#' Method sData<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod sData<-
setGeneric("sData<-", function(object, value) standardGeneric("sData<-"))

#' Method History.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod History
setGeneric("History", function(object) standardGeneric("History"))

#' Method History<-.
#'
#' @name RepSeqExperiment-class
#' @rdname RepSeqExperiment-class
#' @exportMethod History<-
setGeneric("History<-", function(object, value) standardGeneric("History<-"))

#------------------------------------------------------------------
# get methods
#------------------------------------------------------------------

# \code{assay} get count data
# @title The method assay is defined in the class RepSeqExperiment
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]   	
# @return a data.table of assay of VpJ features (row) across samples (columns).
# @name RepSeqExperiment
# @rdname RepSeqExperiment
# @exportMethod assay

#' @rdname RepSeqExperiment-class
#' @aliases assay
#' @aliases assay,RepSeqExperiment-method
setMethod(f = "assay",
    signature = "RepSeqExperiment", 
    definition = function(object) object@assayData
)

# \code{assay<-} set count data
# @title The method assay<- is defined in the class RepSeqExperiment
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]   	
# @exportMethod assay<-
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases assay<-
#' @aliases assay<-,RepSeqExperiment-method
setReplaceMethod(f = "assay",
    signature = "RepSeqExperiment", 
    definition = function(object, i, j, value) {
        set(object@assayData, i, j, value)
        #object@assayData <- value
        object
        }
)

#------------------------------------------------------------------
# meta data
#------------------------------------------------------------------
#' @rdname RepSeqExperiment-class
#' @aliases mData
#' @aliases mData,RepSeqExperiment-method
setMethod(f = "mData",
    signature = "RepSeqExperiment",
    definition = function(object) object@metaData
)


#' @rdname RepSeqExperiment-class
#' @aliases mData<-
#' @aliases mData<-,RepSeqExperiment-method
setReplaceMethod(f = "mData",
    signature = "RepSeqExperiment",
    definition = function(object, value) {
        object@metaData <- value
    object
    }
)

#------------------------------------------------------------------
# sample data
#------------------------------------------------------------------

# \code{sData} get sample information data
# @title The method sData is defined in the class RepSeqExperiment
# @param object an RepSeqExperiment object [\code{\linkS4class{RepSeqExperiment}}]   	
# @return a data frame of sample information, samples are in rows and parameters are in columns.
# @exportMethod sData
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases sData
#' @aliases sData,RepSeqExperiment-method
setMethod(f = "sData",
    signature = "RepSeqExperiment",
    definition = function(object) object@sampleData
)

# \code{sData<-} update sample data 
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}] 
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]   	
# @exportMethod sData<-
# @rdname RepSeqExperiment-class


#' @rdname RepSeqExperiment-class
#' @aliases sData<-
#' @aliases sData<-,RepSeqExperiment-method
setReplaceMethod(
    f = "sData",
    signature = "RepSeqExperiment",
    definition = function(object, value) {
        object@sampleData <- value
    object
    }
)
#------------------------------------------------------------------
# History
#------------------------------------------------------------------
# get history of the object
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}] 
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]   	
# @return a data frame of annotation of clonotype features.
# @exportMethod History
# @rdname RepSeqExperiment


#' @rdname RepSeqExperiment-class
#' @aliases History
#' @aliases History,RepSeqExperiment-method
setMethod(f = "History",
    signature = "RepSeqExperiment",
    definition = function(object) object@History
)

# set history
# @title The method History is defined in the class [\code{\linkS4class{RepSeqExperiment}}] 
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]   	
# @exportMethod History
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases History<-
#' @aliases History<-,RepSeqExperiment-method
setReplaceMethod(f = "History", 
    signature = "RepSeqExperiment", 
    definition = function(object, value) {
        object@History <- value
        object
        }
)

#------------------------------------------------------------------
# display object
#------------------------------------------------------------------
# display the object
# @title The method show is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
# @param object an object of class [\code{\linkS4class{RepSeqExperiment}}]    	
# @exportMethod show
# @rdname RepSeqExperiment

#' @rdname RepSeqExperiment-class
#' @aliases show,RepSeqExperiment-method
setMethod("show", "RepSeqExperiment",
function(object) {
    cts <- assay(object)
    sNames <- unique(cts$lib)
    n <- length(sNames)
    m <- cts[, uniqueN(VpJ)]
    V <- cts[, unique(V)]
    J <- cts[, unique(J)]
    VJ <- cts[, uniqueN(VJ)]
	cat("An object of class \"", class(object), "\"\n", sep="")
	cat("Dimension                  :", m, "clonotypes,", n, "samples\n")
    cat("Number of V genes          :", length(V), "-", V[1:3], "...", V[length(V)], "\n")
	cat("Number of J genes          :", length(J), "-", J[1:3], "...", J[length(J)], "\n")
	cat("Number of V-J genes        :", VJ, "\n")
	#cat("Number of peptide sequences:", length(unique(cts$CDR3aa)), "\n")
	if (n < 4) {
	   cat("Sample names               :", sNames[1:n], "\n")
	   } else {
	       cat("Sample names               :", sNames[1:3], "...", sNames[n],"\n")
	   }
})

#' @rdname RepSeqExperiment-class
#' @aliases names
#' @aliases names,RepSeqExperiment-method
setMethod(f = "names",
    signature(x="RepSeqExperiment"),
    definition = function(x) {
        rownames(sData(x))
    }
)

#' @rdname RepSeqExperiment-class
#' @aliases names
#' @aliases names,RepSeqExperiment-method
setReplaceMethod(f = "names",
    signature(x="RepSeqExperiment", value="ANY"),
    definition = function(x, value) {
        oldnames <- rownames(sData(x))
        rownames(sData(x)) <- value
        snames <- unique(assay(x)[["lib"]])
        for (l in 1:length(snames)) {
            set(assay(x), i=which(assay(x)[["lib"]] == snames[l]), j="lib", value=value[l])
        }
    History(x) <- data.frame(rbind(History(x), 
                data.frame(history=paste(date(), "- updated sample names", paste0(snames, collapse=", "), "from", paste0(snames, collapse=", "), "using names()"))))
    x
    }
)

# set validity
setValidity("RepSeqExperiment", function(object) {
    msg <- NULL
    valid <- TRUE
    if (!identical(unique(assay(object)$lib), rownames(sData(object)))) {
        valid <- FALSE
        msg <- c(msg, "lib column in countData must contain sampleData row names.")
    }
    if (!any(assay(object)$count %% 1 == 0)) {
        valid <- FALSE
        msg <- c(msg, "some count in assay are not integers.")
    }
    if (valid) TRUE else msg
})


#------------------------------------------------------------------
# Subsetting
#------------------------------------------------------------------

# @title The method [ is defined in the class [\code{\linkS4class{RepSeqExperiment}}]
#' Wrapper functions
#'
#' \code{[} Extract parts of RepSeqExperiment.
#' @param i indice(s) of clonotype(s) to extract
#' @param j indice(s) of sample(s) to extract
#' @param drop If \code{TRUE} the result is coerced to the lowest possible dimension
#' @return an object of class [\code{\linkS4class{RepSeqExperiment}}]
#' @rdname RepSeqExperiment-class
#' @aliases [,RepSeqExperiment-method
#' @export
setMethod(
    f = "[",
    signature(x = "RepSeqExperiment", i = "ANY", j = "ANY"),
    definition = function(x, i, j, drop) {
    	if (missing(j)) {
    		j <- 1:nrow(sData(x))
    	}
    	if (!is.character(j)) s <- rownames(sData(x))[j] else s <- j
        cts <- copy(assay(x))
        cts <- cts[lib %in% s, ]
    	out <- new("RepSeqExperiment", 
    				assayData=cts,
    				sampleData=droplevels(sData(x)[j, , drop=FALSE]),
    				metaData=mData(x), 
    				History=data.frame(rbind(History(x),    
    				    data.frame(history=paste0("subet by [ ", length(j)," samples were selected from orignal object RepSeqExperiment."), stringsAsFactors=FALSE)))
	   			)
	   out
})

#' \code{is.RepSeqExperiment} check whether an object is [\code{\linkS4class{RepSeqExperiment}}]
#' @param x an object
# @return TRUE if x is an object of class [\code{\linkS4class{RepSeqExperiment}}].
#' @name is.RepSeqExperiment
#' @rdname RepSeqExperiment-class
#' @export
is.RepSeqExperiment <- function(x) inherits(x, "RepSeqExperiment")
