#' Example location data for spatial genomic data for two data types
#'
#' An example dataset with location for the geographical regions, assume both data types have the same location information for the samples
#'
#' @docType data
#'
#' @usage data(location)
#'
#' @format A data matrix with 34 samples and 2 varialbes:
#' \describe{
#' \item{Column}{The column information for each sample}
#' \item{Row}{The row information for each sample}
#' }
#'
#' @examples
#' data(location)
"location"

#' Example label information for spatial genomic data for two data types
#'
#' An example dataset with pathology subtypes/disease grades for the geographical regions, assume both data types have the same subtype information for the samples
#'
#' @docType data
#'
#' @usage data(label)
#'
#' @format A vector of 34 samples:
#' 
#' @examples
#' data(label)
"label"


#' Example ratio data for spatial genomic data for data type 1
#'
#' An example dataset with the log2ratio compared to healthy controls for the geographical regions for data type 1
#'
#' @docType data
#'
#' @usage data(ratio1)
#'
#' @format A data matrix with 100 genes and 34 samples, rows are gene names and columns are the geographical region/sample names
#'
#' @examples
#' data(ratio1)
"ratio1"

#' Example ratio data for spatial genomic data for data type 2
#'
#' An example dataset with the log2ratio compared to healthy controls for the geographical regions for data type 2
#'
#' @docType data
#'
#' @usage data(ratio2)
#'
#' @format A data matrix with 100 genes and 34 samples, rows are gene names and columns are the geographical region/sample names
#'
#' @examples
#' data(ratio2)
"ratio2"

