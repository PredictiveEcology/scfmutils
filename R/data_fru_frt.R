#' Attributes and fire metrics of the 60 Fire Regime Units (FRUs)
#'
#' Table A1 from Erni et al. (2020).
#'
#' @format
#' A data frame with 60 rows and 10 columns:
#' \describe{
#'   \item{FRU}{Fire Regime Unit polygon id}
#'   \item{Area}{Area (\eqn{Mha})}
#'   \item{Population}{Population (\eqn{no. · 100 ha^{-1}})}
#'   \item{Fire.size}{Fire size (\eqn{ha})}
#'   \item{Frequency}{Fire frequency (\eqn{no. · Mha^{-1} year^{-1}})}
#'   \item{Burn.rate}{Burn rate (\eqn{\% year^{-1}})}
#'   \item{Spring}{Percentage of fires occurring in spring (\eqn{\%})}
#'   \item{Summer}{Percentage of fires occurring in summer (\eqn{\%})}
#'   \item{Lightning}{Percentage of fires caused by lightning (\eqn{\%})}
#'   \item{Human}{Percentage of fires caused by humans (\eqn{\%})}
#' }
#'
#' @note Spring and summer columns provide the distribution, in percent, of area burned during
#' the fire season, i.e. during spring (1 April to 21 June) and summer (22 June to 30 September).
#' Lightning and human columns provide the distribution, in percent, of the number of fires
#' depending on the cause.
#' The median of each fire metric was calculated by year, using fire data ≥ 50 ha for 1970–2016,
#' and these were then averaged to obtain one value per FRU.
#' Bold values indicate the three highest items for each category.
#' `NA`, not available.
#'
#' @source Erni, S., Wang, X., Taylor, S., Boulanger, Y., Swystun, T., Flannigan, M.,
#' Parisien, M.-A.. 2020. Developing a two-level fire regime zonation system for Canada.
#' Canadian Journal of Forest Research. 50(3): 259-273. <https://doi.org/10.1139/cjfr-2019-0191>
#'
"fru_attr"

#' Attributes and fire metrics of the 15 Fire Regime Units (FRTs)
#'
#' Table A2 from Erni et al. (2020).
#'
#' @format
#' A data frame with 15 rows and 10 columns:
#' \describe{
#'   \item{FRT}{Fire Regime Type polygon id}
#'   \item{Area}{Area (\eqn{Mha})}
#'   \item{Population}{Population (\eqn{no. · 100 ha^{-1}})}
#'   \item{Fire.size}{Fire size (\eqn{ha})}
#'   \item{Frequency}{Fire frequency (\eqn{no. · Mha^{-1} year^{-1}})}
#'   \item{Burn.rate}{Burn rate (\eqn{\% year^{-1}})}
#'   \item{Spring}{Percentage of fires occurring in spring (\eqn{\%})}
#'   \item{Summer}{Percentage of fires occurring in summer (\eqn{\%})}
#'   \item{Lightning}{Percentage of fires caused by lightning (\eqn{\%})}
#'   \item{Human}{Percentage of fires caused by humans (\eqn{\%})}
#' }
#'
#' @note Spring and summer columns provide the distribution, in percent, of area burned during
#' the fire season, i.e. during spring (1 April to 21 June) and summer (22 June to 30 September).
#' Lightning and human columns provide the distribution, in percent, of the number of fires
#' depending on the cause.
#' The median of each fire metric was calculated by year, using fire data ≥ 50 ha for 1970–2016,
#' and these then averaged to obtain one value per FRT.
#' Bold values indicate the three highest items for each category.
#'
#' @source Erni, S., Wang, X., Taylor, S., Boulanger, Y., Swystun, T., Flannigan, M.,
#' Parisien, M.-A.. 2020. Developing a two-level fire regime zonation system for Canada.
#' Canadian Journal of Forest Research. 50(3): 259-273. <https://doi.org/10.1139/cjfr-2019-0191>
#'
"frt_attr"
