#' Retrieve CDEC CSV
#'
#' @param sta 3-character alphanumeric code for CDEC station
#' @param sen Numeric code for CDEC sensor
#' @param dur 1-character code for duration type: `D`=daily, `H`=hourly, `E`=event (15 minute)
#' @param start_date Start date for query expressed as R date object e.g. lubridate::ymd("1900-10-01")
#' @param end_date End date for query expressed as R date object e.g. lubridate::ymd("2023-09-30")
#' @param dir Directory to store outputs in
#'
#' @returns The path to the downloaded CDEC file in `csv.gz` (compressed CSV) format
#' @md
#'
#' @examples
#' @export
cdec_csv_retrieve <- function(sta=character(), sen=character(), dur=character(),
                              start_date=lubridate::ymd("1900-10-01"), end_date=lubridate::ymd("2023-09-30"),
                              dir="temp") {
  name <- str_to_upper(paste0(sta, "_", sen, "_", dur))
  dir.create(dir, recursive=T)
  filename <- file.path(dir, paste0(name,".csv.gz"))
  if(!file.exists(filename)){
    message(paste0("downloading to ",filename))
    dir.create(dir, recursive = TRUE)
    data_raw <-
      httr::GET(
        url="https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet",
        query=list(Stations=sta, SensorNums=sen, dur_code=dur,
                   Start=format(start_date,"%Y-%m-%d"), End=format(end_date,"%Y-%m-%d"))) |>
      httr::content("raw")
    gzf <- gzcon(file(filename, "wb"))
    data_raw |> writeBin(gzf)
    close(gzf)
  } else {
    message(paste0(filename, " exists, loading"))
  }
  return(filename)
}

#' Parse CDEC CSV
#'
#' @param filename Path to a downloaded CDEC file in the format output by the CDEC `CSVDataServlet` service, such as retrieved by the `retrieve_cdec_csv()` function. May be compressed (`.csv.gz`) or uncompressed (`.csv`).
#'
#' @returns A `tbl_df` data frame with cleaned tabular data
#' @md
#'
#' @examples
#' @export
cdec_csv_parse <- function(filename) {
  read_csv(filename,
           col_select = c(date_time = "DATE TIME", value = "VALUE"),
           col_types = list(col_datetime(), col_double())) |>
    janitor::clean_names() |>
    mutate(water_year = if_else(month(date_time)>=10, year(date_time)+1, year(date_time)))
}

#' Calculate Days Inundated
#'
#' @param q_gauge A vector of flow values, typically a time series
#' @param q_crit A numeric value against which to calculate the exceedance
#' @param stat Either "max" to calculate maximum days inundated (default) or "sum" to calculate total days inundated
#'
#' @returns A single numeric value for the specified statistic.
#' @md
#'
#' @examples
#' @export
duration_calc_days_inundated <- function(q_gauge, q_crit, stat="max") {

  n_days <- length(q_gauge)

  exceedence_intervals <-
    tibble(exceeds_q = q_gauge > q_crit,
           exceedence_event = with(rle(exceeds_q), rep(seq_along(lengths), lengths))) |>
    group_by(exceedence_event) |>
    mutate(cumulative_exceedence = if_else(exceeds_q, seq_along(exceedence_event), 0)) |>
    filter(exceeds_q) |>
    summarize(duration_days = max(cumulative_exceedence))

  if(stat=="sum"){
    return(sum(exceedence_intervals$duration_days))
  } else if(stat=="max") {
    return(pmax(max(exceedence_intervals$duration_days),0))
  }
}

#' Apply Duration HSI to Flow-to-Suitable-Area Curve
#'
#' @param fsa A `data.frame` or `tbl_df` containing a flow-to-suitable-area curve
#' @param drc A `data.frame` or `tbl_df` containing a flow-to-duration-HSI curve
#' @param fsa_q The unquoted name of the flow column in `fsa` (defaults to `q`)
#' @param fsa_wua The unquoted name of the suitable area column in `fsa` (defaults to `wua`)
#' @param drc_q The unquoted name of the flow column in `drc` (defaults to `q`)
#' @param drc_dhsi The unquoted name of the duration suitability factor column in `drc` (defaults to `dhsi`)
#'
#' @return A `tbl_df` data frame with one row per flow (`q`) containing columns for the original suitable area (`wua`) and the new duration-weighted suitable area (`durwua`) along with interim calculation columns.
#' @md
#'
#' @examples
#' @export
duration_apply_dhsi_to_fsa_curve <- function(fsa, drc,
                                         fsa_q=q, fsa_wua=wua,
                                         drc_q=q, drc_dhsi=dhsi) {

  fun <- function(fsa, drc,
                  fsa_q={{fsa_q}}, fsa_wua={{fsa_wua}},
                  drc_q={{drc_q}}, drc_dhsi={{drc_dhsi}}) {

    full_join(fsa |> transmute(q={{fsa_q}}, wua={{fsa_wua}}),
              drc |> transmute(q={{drc_q}}, dhsi={{drc_dhsi}}),
              by = join_by(q)) |>
      # fill in any gaps with linear interpolation
      arrange(q) |>
      mutate(wua = zoo::na.approx(wua, q, na.rm=F),
             dhsi = zoo::na.approx(dhsi, q, na.rm=F)) |>
      drop_na() |>
      # start by computing the marginal WUA added at each step
      mutate(marginal_wua = if_else(row_number()==1, wua, wua - lag(wua))) |>
      # apply each flow's DHSI to whatever area is added or removed at that flow
      ### mutate(marginal_durwua = pmin(dhsi * marginal_wua, marginal_wua),
      ###       # pmin is a simplification to account for loss; this will end up overestimating the magnitude of the marginal drop in habitat
      ###       durwua = pmax(cumsum(marginal_durwua), 0))
      mutate(wua_added = if_else(marginal_wua > 0, marginal_wua, 0),
             marginal_durwua_added = dhsi * wua_added,
             wua_removed = (-1) * if_else(marginal_wua < 0, marginal_wua, 0),
             marginal_durwua_removed = cummean(dhsi) * wua_removed,
             ###        marginal_durwua_removed = cummean.na.rm(if_else(dhsi > 0, dhsi, NA)) * wua_removed,
             marginal_durwua = marginal_durwua_added - marginal_durwua_removed,
             durwua = pmin(pmax(cumsum(marginal_durwua), 0), wua))
  }

  fun <- possibly(fun, otherwise=NA) # apply error handling

  return(fun(fsa, drc,
             fsa_q={{fsa_q}}, fsa_wua={{fsa_wua}},
             drc_q={{drc_q}}, drc_dhsi={{drc_dhsi}}))
}

