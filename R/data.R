#' New York Auto Insurance Complaint Rankings
#'
#' Automobile insurance company complaint data from the New York State
#' Department of Financial Services, covering filing years 2009--2020. Each row
#' represents one insurer in one filing year.
#'
#' @format A data frame with 1942 rows and 4 variables:
#' \describe{
#'   \item{year}{Filing year (2009--2020).}
#'   \item{upheld}{Number of upheld complaints against the insurer.}
#'   \item{total}{Total complaints filed against the insurer.}
#'   \item{premiums}{Premiums written (millions of USD).}
#' }
#'
#' @source New York State Department of Financial Services,
#'   Automobile Insurance Company Complaint Rankings, Beginning 2009.
"ny_complaints"

#' Doctor Visits
#'
#' Number of doctor visits for 1812 individuals from the German
#' Socio-Economic Panel (GSOEP). Originally distributed in the \pkg{zic}
#' package.
#'
#' @format A data frame with 1812 rows and 1 variable:
#' \describe{
#'   \item{visits}{Number of doctor visits (non-negative integer).}
#' }
#'
#' @source Winkelmann, R. (2004). Health care reform and the number of doctor
#'   visits -- an econometric analysis. \emph{Journal of Applied Econometrics},
#'   19(4), 455--472.
"docvisits"

#' Gaming and Betting Offenses in New South Wales
#'
#' Counts of gaming and betting offenses recorded at 342 locations in
#' New South Wales, Australia.
#'
#' @format A data frame with 342 rows and 1 variable:
#' \describe{
#'   \item{offenses}{Number of gaming and betting offenses (non-negative integer).}
#' }
#'
#' @source Australian Bureau of Statistics, New South Wales crime data.
"nsw_offenses"
