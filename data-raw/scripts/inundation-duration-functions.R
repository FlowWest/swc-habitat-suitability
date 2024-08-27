# INUNDATION DURATION ##########################################################

# STEP 1: CDEC Retrieval -------------------------------------------------------
# Retrieve water year gauge data from CDEC
# functions from interoperable flows

retrieve_cdec_csv <- function(sta=character(), sen=character(), dur=character(),
                              start_date=ymd("1900-10-01"), end_date=ymd("2023-09-30"),
                              dir="temp") {
  name <- str_to_upper(paste0(sta, "_", sen, "_", dur))
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

parse_cdec_csv <- function(filename) {
  read_csv(filename,
           col_select = c(date_time = "DATE TIME", value = "VALUE"),
           col_types = list(col_datetime(), col_double())) |>
    janitor::clean_names() |>
    mutate(water_year = if_else(month(date_time)>=10, year(date_time)+1, year(date_time)))
}

# TODO: Function to post-process CDEC gauge data to get it ready for next step

# STEP 2: Days Inundated By Flow -----------------------------------------------
# Calculate number of days inundated per water year at each flow in a list
# The list of flows is the flows of a flow-to-suitable-area curve

calculate_days_inundated <- function(q_gauge, q_crit, stat="max") {

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

# TODO: Function to turn this duration result into a duration HSI curve

# STEP 3: Apply Duration HSI Curve to Flow-to-Suitable-Area Curve --------------

cummean.na.rm <- function(x) { # https://stackoverflow.com/a/73215592
  idx <- cumsum(!is.na(x))
  x_filtered <- x[!is.na(x)]
  cummean(x_filtered)[idx]
}

apply_durations_to_fsa_curve <- function(fsa, drc, fsa_q=q, fsa_wua=wua, drc_q=q, drc_dhsi=dhsi, monotonic=F) {
  #tryCatch(
  #  expr = {
      m <- if (monotonic) 0 else 1 # forces monotonic if user has selected it so

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
               marginal_durwua_removed = cummean(dhsi) * wua_removed * m,
               ###        marginal_durwua_removed = cummean.na.rm(if_else(dhsi > 0, dhsi, NA)) * wua_removed,
               marginal_durwua = marginal_durwua_added - marginal_durwua_removed, # to modify so that it doesn't drop negative before habitat has even started being added, use:
               #marginal_durwua = if_else(cumsum(marginal_durwua_added)>0, marginal_durwua_added - marginal_durwua_removed, 0),
               durwua_raw = pmin(pmax(cumsum(marginal_durwua), 0), wua),
               marginal_durwua_bounded = if_else((cumsum(marginal_durwua) > 0) | (marginal_durwua > 0), marginal_durwua, 0),
               #durwua = pmin(pmax(cumsum(marginal_durwua), 0), wua)
               durwua = pmin(pmax(cumsum(marginal_durwua_bounded), 0), wua)
               )
  #  },
  #error = function(e) {
  #  message(e)
  #  return(NA)
  #})
}
apply_durations_to_fsa_curve <- possibly(apply_durations_to_fsa_curve, otherwise=NA)

# Baseflow Elimination ---------------------------------------------------------

eliminate_baseflow <- function(fsa, baseflow_threshold_flow, fsa_q=q, fsa_wua=wua, ...) {

  min_q <- 0
  max_q <- fsa |> pull({{fsa_q}}) |> max()

  baseflow_removal_drc <-
    fsa |>
    complete(flow_cfs = baseflow_threshold_flow) |>
    arrange(flow_cfs) |>
    transmute(q = flow_cfs,
              dhsi = if_else(flow_cfs <= baseflow_threshold_flow, 0, 1))

  apply_durations_to_fsa_curve(fsa=fsa, drc=baseflow_removal_drc,
                               fsa_q={{fsa_q}}, fsa_wua={{fsa_wua}}, ...) |>
   transmute({{fsa_q}} := q,
             {{fsa_wua}} := durwua)

}

