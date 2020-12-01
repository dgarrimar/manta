#' Simulated Measurements of Five Disease Biomarkers
#'
#' A simulated dataset containing the levels of 5 biomarkers, 
#' measured in 100 individuals, with different scales. 
#' Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 5 numerical variables:
#' 
#' \describe{
#'   \item{biomarker1}{levels of biomarker1}
#'   \item{biomarker2}{levels of biomarker2}
#'   \item{...}{}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name biomarkers
#' @usage data(biomarkers)
#' @author Diego Garrido-Martín
#' 
NULL

#' Simulated Metadata for 100 Patients
#'
#' A simulated dataset containing the age, gender and disease status of 100
#' individuals. Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 3 variables:
#' 
#' \describe{
#'   \item{age}{Age of the patient (numerical)}
#'   \item{gender}{Gender of the patient (factor with levels: "male" and 
#'   "female")}
#'   \item{status}{Disease status of the patient (ordered factor with levels:
#'   "healthy", "mild" and "severe")}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name patients
#' @usage data(patients)
#' @author Diego Garrido-Martín
#' 
NULL