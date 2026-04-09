#' Daily Shallow Cumulus and Atmospheric Boundary Layer Dataset from ARF
#'
#' The `arf_daily` dataset contains daily-aggregated atmospheric, surface,
#' and boundary layer variables derived from ceilometer and meteorological
#' observations. The dataset is designed for analysis of boundary layer
#' structure, cloud properties, and surface-atmosphere interactions on a
#' daily temporal scale. For all variables below, 'daytime' referes to
#' 8am-8pm local time.
#'
#' @format A data frame with \code{n} rows (days) and \code{p} variables:
#' \describe{
#'   \item{date}{Date. Calendar date of observation, in unix time, format 'YYYY-MM-DD'.}
#'
#'   \item{pblh_RZ_q}{Numeric. Planetary boundary layer height estimated using humidity-based Richardson number method (m).}
#'   \item{pblh_RZ_thta}{Numeric. Boundary layer height estimated using potential temperature Richardson number method (m).}
#'   \item{pblh_parcel}{Numeric. Boundary layer height estimated using parcel method (m).}
#'   \item{pblh_grad_thta}{Numeric. Boundary layer height derived from potential temperature gradient method (m).}
#'   \item{pblh_gradq}{Numeric. Boundary layer height derived from humidity gradient method (m).}
#'
#'   \item{median_PBLH}{Numeric. Median boundary layer height across all 5 methods (m).}
#'   \item{avg_PBLH}{Numeric. Mean boundary layer height across all 5 methods (m).}
#'
#'   \item{WD_1_2_1_daytime_mean}{Numeric. Mean circular wind direction during daytime (degrees).}
#'   \item{WS_1_2_1_daytime_mean}{Numeric. Mean wind speed during daytime (m/s).}
#'
#'   \item{T_CANOPY_daytime_mean}{Numeric. Mean canopy temperature during daytime (°C).}
#'   \item{PA_1_2_1_daytime_mean}{Numeric. Mean atmospheric pressure during daytime, retrieved from lowest tower level (2m) (kPa).}
#'   \item{TA_1_2_1_daytime_mean}{Numeric. Mean air temperature during daytime, retrieved from lowest tower level (2m) (°C).}
#'   \item{TA_1_7_1_daytime_mean}{Numeric. Mean air temperature during daytime, retrieved from highest tower level (28m) (°C).}
#'
#'   \item{RH_1_2_1_daytime_mean}{Numeric. Mean relative humidity at level 1 during daytime, retrieved from lowest tower level (2m) (\%).}
#'   \item{RH_1_7_1_daytime_mean}{Numeric. Mean relative humidity at level 7 during daytime, retrieved from highest tower level (28m) (\%).}
#'
#'   \item{SWC_4_1_1_daytime_mean}{Numeric. Mean soil water content at location 2, shallowest sensor (m³/m³).}
#'   \item{SWC_4_2_1_daytime_mean}{Numeric. Mean soil water content at location 2, instermediate sensor (m³/m³).}
#'   \item{SWC_4_3_1_daytime_mean}{Numeric. Mean soil water content at location 2, deepest sensor (m³/m³).}
#'   \item{SWC_3_1_1_daytime_mean}{Numeric. Mean soil water content at location 1, shallowest sensor (m³/m³).}
#'
#'   \item{TS_PI_F_1_1_1_daytime_mean}{Numeric. Mean soil temperature (°C).}
#'
#'   \item{daytime_hours}{Numeric. Total number of daylight hours for the day, calculated using "astral" python package (hours).}
#'   \item{year}{Integer. Year of observation.}
#'   \item{doy}{Integer. Day of year (1–366).}
#'
#'   \item{LE_f_daytime_sum}{Numeric. Daytime sum of latent heat flux (W/m²/day).}
#'   \item{H_f_daytime_sum}{Numeric. Daytime sum of sensible heat flux (W/m²/day).}
#'
#'   \item{SW_IN_daytime_sum}{Numeric. Incoming shortwave radiation summed over daytime hours (W/m²/day).}
#'
#'   \item{G_1_1_1_daytime_sum}{Numeric. Ground heat flux (sensor 1) summed over daytime (W/m²/day).}
#'   \item{G_2_1_1_daytime_sum}{Numeric. Ground heat flux (sensor 2) summed over daytime (W/m²/day).}
#'
#'   \item{blc_flag}{Integer/Logical. Indicator for boundary layer coupling state ([0 = no BLC], [1 = Clean BLC], [2 = Messy/not ideal BLC]).}
#'   \item{Conditions}{Integer/Logical. Integer value indicating daily boundary layer conditions, derived from ceilomter and satellite imagery ([0 = clear sky], [1 =  Morning clear sky to afternoon shallow cumulus transition], [2 = stratiform/overcast], [3 = Morning overcast to afternoon shallow cumulus transition], [4 = synoptic weather (rain/snow/fog)], [5 = Morning synoptic weather to afternoon shallow cumulus transition]).}
#'   \item{Smoke}{Binary. Binary column indicating presence of smoke or aerosol in atmospheric boundary layer. 1 = smoke, 0 = no smoke.}
#'
#'   \item{LCL_height_km}{Numeric. Lifting condensation level, derived from Caribou soundings (km).}
#'   \item{LTS_K}{Numeric. Lower tropospheric stability (i.e., change in potential temperature between 700 hPa level and surface), derived from caribou soundings (K).}
#'   \item{LapseRate_850_500_CperKm}{Numeric. Lapse rate between 850 and 500 hPa, derived from Caribou soundings (°C/km).}
#'   \item{PW_mm}{Numeric. Precipitable water (integrated depth of condensed water in vertical column) (mm).}
#'
#'   \item{mean_cf_2_6}{Numeric. Mean cloud fraction between 02–06pm local time, decimal between 0-1.}
#'   \item{std_cf_2_6}{Numeric. Standard deviation of cloud fraction between 02–06pm local time.}
#'   \item{cv_cf_2_6}{Numeric. Coefficient of variation of cloud fraction between 02–06pm local time.}
#'
#'   \item{mean_cbh_2_6}{Numeric. Mean cloud base height between 02–06pm local time, derived from ceilometer (m).}
#'   \item{std_cbh_2_6}{Numeric. Standard deviation of cloud base height 02–06pm local time, derived from ceilometer (m).}
#'
#'   \item{morning_cbh}{Numeric. Morning (8am-11am LST averaged) cloud base height (m).}
#'   \item{evening_cbh}{Numeric. Evening (4pm-7pm LST averaged) cloud base height (m).}
#'   \item{cbh_rise}{Numeric. Change in cloud base height from morning to evening (m).}
#'   \item{cbh_rise_filtered}{Numeric. Filtered change in CBH after quality control, eliminating all non-BLC days using blc_flag (m).}
#'
#'   \item{mean_LCL_4_7}{Numeric. Mean LCL between 04–07 local time, calculated using metpy, from local tower observations (vs LCL_height_km from Caribou soundings) (m).}
#' }
#'
#' @details
#' This dataset integrates ceilometer-derived cloud properties with surface
#' meteorological and flux measurements from a forest site near Fredericton, NB,
#' and supplements this data with nearby weather balloon soundings in Caribou, MA.
#' Boundary layer height (PBLH) is
#' estimated using multiple methods, allowing intercomparison and uncertainty
#' analysis. Cloud metrics are computed over specific time windows to capture
#' peak boundary layer cloud times, and diurnal variability.
#'
#' Daytime averages and sums are computed using local solar time constraints,
#' as we are not interested in nighttime data.
#'
#' @source
#' Derived from ceilometer observations and surface meteorological data.
#'
#' @usage data(arf_daily)
#'
#' @examples
#' data(arf_daily)
#' str(arf_daily)
#' summary(arf_daily)
"arf_daily"
