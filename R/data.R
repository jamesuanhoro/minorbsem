#' Holzinger-Swineford data
#'
#' Reduced form data retrieved from lavaan
#'
#' @format ## `HS`
#' A data frame with 301 observations and 15 variables:
#' \describe{
#'   \item{id}{Identifier}
#'   \item{sex}{Gender}
#'   \item{ageyr}{Age, year part}
#'   \item{agemo}{Age, month part}
#'   \item{school}{School (Pasteur or Grant-White)}
#'   \item{grade}{Grade}
#'   \item{x1}{Visual perception}
#'   \item{x2}{Cubes}
#'   \item{x3}{Lozenges}
#'   \item{x4}{Paragraph comprehension}
#'   \item{x5}{Sentence completion}
#'   \item{x6}{Word meaning}
#'   \item{x7}{Speeded addition}
#'   \item{x8}{Speeded counting of dots}
#'   \item{x9}{Speeded discrimination straight and curved capitals}
#' }
#' @source <https://cran.r-project.org/web/packages/lavaan/index.html>
"HS"

#' Industrialization And Political Democracy Dataset
#'
#' Dataset in Bollen (1989), also retrieved from lavaan
#'
#' @format ## `PD`
#' A data frame with 75 observations and 11 variables:
#' \describe{
#'   \item{y1}{Expert ratings of the freedom of the press in 1960}
#'   \item{y2}{The freedom of political opposition in 1960}
#'   \item{y3}{The fairness of elections in 1960}
#'   \item{y4}{The effectiveness of the elected legislature in 1960}
#'   \item{y5}{Expert ratings of the freedom of the press in 1965}
#'   \item{y6}{The freedom of political opposition in 1965}
#'   \item{y7}{The fairness of elections in 1965}
#'   \item{y8}{The effectiveness of the elected legislature in 1965}
#'   \item{x1}{The gross national product (GNP) per capita in 1960}
#'   \item{x2}{The inanimate energy consumption per capita in 1960}
#'   \item{x3}{The percentage of the labor force in industry in 1960}
#' }
#' @source <https://cran.r-project.org/web/packages/lavaan/index.html>
"PD"

#' Studies on the Hospital Anxiety and Depression Scale.
#'
#' The data set includes 28 studies on 14 items measuring the Hospital Anxiety
#' and Depression Scale (HADS) Reported by
#' \insertCite{norton_hospital_2013;textual}{minorbsem}, with data and
#' documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{minorbsem}.
#'
#' @format ## `Norton13`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{
#'      A list of 28 studies of correlation matrices.
#'      The variables are 14 items (x1 to x14) measuring HADS}
#'   \item{n}{A vector of sample sizes}
#'   \item{population}{A vector of the population of the data}
#'   \item{group}{
#'      A vector of classification into patients vs. non-patients
#'      based on population}
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"Norton13"

#' Another dataset from the metaSEM package.
#'
#' The data set includes 11 studies on 9 items measuring work-related attitudes
#' \insertCite{noauthor_international_1992}{minorbsem}.
#' Data and documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{minorbsem}.
#'
#' @format ## `issp89`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{A list of 11 studies of covariance matrices}
#'   \item{n}{A vector of sample sizes}
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"issp89"
