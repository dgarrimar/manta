#' Biomarkers
#'
#' A simulated dataset containing the levels of 5 biomarkers, 
#' measured in 100 individuals, with different scales. 
#' Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 5 numerical variables:
#' \describe{
#'   \item{biomarker1}{levels of biomarker1}
#'   \item{biomarker2}{levels of biomarker2}
#'   ...
#' }
#' @author Diego Garrido-Martín
#' 
"biomarkers"

#' Patients
#'
#' A simulated dataset containing the gender, age and disease status of 100
#' individuals. Missing observations appear as \code{NA}.
#'
#' @format A matrix with 100 rows and 3 variables:
#' \describe{
#'   \item{gender}{Gender of the patient (factor with levels: \code{male} and 
#'   \code{female})}
#'   \item{age}{Age of the patient (numerical)}
#'   \item{disease}{Disease status of the patient (ordered factor with levels
#'   \code{healthy}, \code{"mild"}, \code{"severe"})}
#' }
#' @author Diego Garrido-Martín
#' 
"patients"
