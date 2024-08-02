function(input, output, session){

  # REACTIVE FLOW FILTER -------------------------------------------------------

  message(paste(names(predictions), collapse = ", "))

  active_predictions <- reactive({
    # predictions_table |>
    #   filter(flow_cfs == input$active_flow) |>
    #   mutate(wua_per_lf = !!sym(input$wua_var))
    predictions |>
      filter(flow_cfs == input$active_flow) |>
      filter(habitat == input$habitat_type) |>
      glimpse() #|>
    # TODO could also filter for selected model and bfc removed/not removed
    #message(nrow(result), " result rows")
    #return(result)
  })

  active_geom <- reactive({
    geom |>
      inner_join(active_predictions(),
                 by = join_by(comid)) |>
      mutate(wua_per_lf = !!sym(input$wua_var)) |>
      filter(!is.na(wua_per_lf)) |>
      mutate(object_label = paste0("<strong>", comid, " ", gnis_name, "</strong>", "<br />",
                                   "Scale-Dependent Model: ",
                                   round(wua_per_lf_pred_SD,2), " ft2/ft", "<br />",
                                   "Scale-Dependent Model (post-model baseflow removal): ",
                                   round(wua_per_lf_pred_SD_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model: ",
                                   round(wua_per_lf_pred_SN,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model (post-model baseflow removal): ",
                                   round(wua_per_lf_pred_SN_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Actual: ", round(wua_per_lf_actual,2), " ft2/ft")) |>
      #filter(watershed_name == active_watershed_name()) |>
      glimpse()
  })

  # SELECTED COMID OBSERVER ----------------------------------------------------

  selected_point <- reactiveValues(object_id = NULL,
                                   comid = NULL,
                                   lng = NULL,
                                   lat = NULL)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 6) == "comid_") {
        selected_point$object_id <- input$main_map_shape_click$id
        selected_point$lng <- input$main_map_shape_click$lng
        selected_point$lat <- input$main_map_shape_click$lat
        selected_point$comid <- as.numeric(str_replace(selected_point$object_id, "comid_", ""))
      }
    }
  })

  active_pred_table <- reactive({
    if (!is.null(selected_point$object_id)) {
      active_predictions() |>
        filter(comid == selected_point$comid) |>
        select(where(is.numeric)) |>
        pivot_longer(everything()) |>
        mutate(across(value, function(x) round(x, 2)))
    }
  })

  output$pred_table <- DT::renderDT({
    if (!is.null(selected_point$object_id)) {
      DT::datatable(active_pred_table(),
                    options = list(paging = F, searching = F))
    }
  })

  active_attr_table <- reactive({
   #if (!is.null(selected_point$object_id) & (length(selected_point$comid)>0)) {
     if (!is.null(selected_point$object_id)) {
     result <- attr |>
       filter(comid == selected_point$comid) |>
       select(where(is.numeric)) |>
       pivot_longer(everything())
       #gc()
     return(result)
   }
  })

  output$attr_table <- DT::renderDT({
      if (!is.null(selected_point$object_id)) {
        DT::datatable(active_attr_table(),
                      options = list(paging = F, searching = F))
      }
    })

  output$fsa_plot <- renderPlot({
    if (!is.null(selected_point$object_id)) { #& (length(selected_point$comid)>0)) {
    predictions |>
      filter(comid == selected_point$comid) |>
      filter(habitat == input$habitat_type) |>
      ggplot(aes(x = flow_cfs)) +
      geom_line(aes(y = wua_per_lf_actual, color="Actual", linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal")), linewidth=2) +
      geom_line(aes(y = wua_per_lf_pred_SD, color="Scale-Dependent", linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
      geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
      geom_line(aes(y = wua_per_lf_pred_SN, color="Scale-Normalized", linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
      geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Analysis")) +
      #geom_hline(aes(yintercept = chan_width_ft)) + #, linetype="Channel Width (ft)")) +
      #geom_vline(aes(xintercept = baseflow_cfs)) +
      #geom_text(aes(x = 1, y = chan_width_ft, label = chan_width_ft)) +
      scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
      scale_y_continuous(limits = c(0, NA)) +
      #scale_y_continuous(trans = ihs, labels = scales::label_comma(), limits = c(0, NA)) +
      theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical") +
      xlab("Flow (cfs)") + ylab("WUA (ft2) per linear ft") +
      scale_color_discrete(name = "Model Type") +
      scale_linetype_manual(name = "Baseflow Method",
                            values = c("Prior BFC Removal" = "solid",
                                       "No BFC Removal" = "solid",
                                       "Duration Analysis" = "dashed",
                                       "Post-Model BFC Removal" = "dotted"))
  } else {
      ggplot()
    }})

  # LEAFLET MAP FUNCTIONS ------------------------------------------------------

  make_leaflet <- function(bbox=c(xmin=-122.3, ymin=38.5, xmax=-121.3, ymax=39.7)) {
    m <- leaflet::leaflet() |>
      leaflet::addMapPane("Basemap", zIndex = 400) |>
      leaflet::addMapPane("Watersheds", zIndex = 440) |>
      leaflet::addMapPane("Flowlines", zIndex = 470) |>
      leaflet::addMapPane("AOI", zIndex = 480) |>
      leaflet::addMapPane("Reference", zIndex = 490) |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 1,
                                                       maxNativeZoom = 13,
                                                       maxZoom = 13,
                                                       pane = "Basemap"
                        )) |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Reference_Overlay/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 1,
                                                       #minZoom = 10,
                                                       pane = "Reference"
                        )) |>
      leaflet::addTiles(urlTemplate = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
                        options = leaflet::tileOptions(noWrap = TRUE,
                                                       opacity = 0.5,
                                                       minZoom = 14,
                                                       maxNativeZoom = 23,
                                                       pane = "Basemap"
                        )) |>
      leaflet::addTiles(urlTemplate = "", attribution = 'Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, IGN, IGP, UPR-EGP, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, Esri, and the GIS User Community') |>
      leaflet::fitBounds(lng1 = bbox[["xmin"]],
                         lat1 = bbox[["ymin"]],
                         lng2 = bbox[["xmax"]],
                         lat2 = bbox[["ymax"]])
    return(m)
  }

  layer_flowlines <- function(m, show = TRUE) {
    if(show) {
      #pal_limits <- c(0, 100) #,
      pal_limits <- c(min(active_geom()$wua_per_lf), max(active_geom()$wua_per_lf))
      message("plotting ", nrow(active_geom()), " flowlines")
      #m |> leafgl::addGlPolylines(data = active_geom(),layerId = ~object_id)
      m |> leaflet::addPolylines(data = active_geom(), # TODO switch to leafgl::addGlPolylines
                                 layerId = ~object_id,
                                 label = ~lapply(object_label, htmltools::HTML),
                                 color = ~pal(wua_per_lf), #pal(wua_per_lf),
                                 opacity = 1,
                                 weight = 2,
                                 options = leaflet::pathOptions(pane = "Flowlines"),
                                 group = "flowlines",
                                 highlightOptions = leaflet::highlightOptions(color = "#8B0000",
                                                                              weight = 3,
                                                                              bringToFront = TRUE)
      ) |>
        leaflet::addLegend(position = "bottomright",
                           colors = rev(pal(seq(pal_limits[[1]], pal_limits[[2]], (pal_limits[[2]]-pal_limits[[1]])/5))),
                           labels = c(round(pal_limits[[2]],2), rep("", 5-1), round(pal_limits[[1]],2)),
                           title = "WUA (ft2) per linear ft",
                           layerId = "clegend")
    } else {
      m |> leaflet::removeShape(active_geom()$comid) |> leaflet::removeControl("clegend")
      # TODO switch to leaflet::removeGlPolylines
    }
  }

  # LEAFLET RENDER -------------------------------------------------------------

  output$main_map <- renderLeaflet({
    shinyjs::showElement(id = 'loading_action')
    make_leaflet() |>
      leaflet::addPolygons(data = watersheds,
                           stroke = T,
                           weight = 0,
                           color = "#FFFFFF",
                           opacity = 0,
                           fill = T,
                           fillColor = "#FFFFFF",
                           fillOpacity = 0,
                           layerId = ~watershed_id,
                           label = ~lapply(watershed_label, htmltools::HTML),
                           highlightOptions = highlightOptions(stroke = T,
                                                               weight = 1,
                                                               color = "white",
                                                               opacity = 1,
                                                               fill = T,
                                                               fillColor = "red",
                                                               fillOpacity = 0.5,
                                                               bringToFront = F)
      ) |>
      layer_flowlines() #|>
    #leaflet::addLayersControl(baseGroups = c("watersheds", "none"),
    #                          overlayGroups = c("flowlines", "aoi"),
    #                          position = "bottomleft",
    #                          options = leaflet::layersControlOptions(collapsed = FALSE))
  })


# SELECTED WATERSHED OBSERVER ----------------------------------------------------

selected_watershed <- reactiveValues(object_id = NA,
                                     lng = NA,
                                     lat = NA)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 10) == "watershed_") {
        selected_watershed$object_id <- input$main_map_shape_click$id
        selected_watershed$lng <- input$main_map_shape_click$lng
        selected_watershed$lat <- input$main_map_shape_click$lat
      }
    }
  })

  active_watershed_geom <- reactive({
    watersheds |>
      filter(watershed_id == selected_watershed$object_id) |>
      mutate(object_id = "active_watershed")
  })

  active_watershed_bbox <- reactive({
    st_bbox(active_watershed_geom())
  })

  active_watershed_name <- reactive({
    if(!is.na(selected_watershed$object_id)) {
      watersheds$watershed_level_3[[which(watersheds$watershed_id==selected_watershed$object_id)]]
    } else {
      NA
    }
  })

  observe({

    message(selected_watershed$object_id)

    proxy <- leaflet::leafletProxy("main_map") |>
      leaflet::removeShape("active_watershed")

    if (nrow(active_watershed_geom()) > 0) {
      proxy |>
        leaflet::addPolygons(data = active_watershed_geom(),
                             stroke = T,
                             weight = 2,
                             color = "red",
                             opacity = 1,
                             fill = T,
                             fillColor = "red",
                             fillOpacity = 0.25,
                             layerId = "active_watershed",
                             label = ~lapply(watershed_label, htmltools::HTML))

      proxy |>
        leaflet::flyToBounds(lng1 = active_watershed_bbox()$xmin,
                             lat1 = active_watershed_bbox()$ymin,
                             lng2 = active_watershed_bbox()$xmax,
                             lat2 = active_watershed_bbox()$ymax)
    }
  })

  ### DURATION ANALYSIS --------------------------------------------------------

  streamgages <- get_data(streamgage_attr, package = "habistat") |>
    transmute(watershed, station_id,
              station_label = glue::glue("{str_to_upper(station_id)}: {str_to_upper(name)} ({min_wy}-{max_wy})")) |>
    nest(.by = watershed) |>
    mutate(data = map(data, function(x) deframe(x) |> as.list())) |>
    deframe()

  streamgage_drc <- reactive({

    active_comid_gradient_class <- (attr |>
                                      filter(comid == selected_point$comid) |>
                                      pull(hqt_gradient_class))[[1]]

    message("streamgage_drc for ", input$streamgage_id,
            " ", input$selected_run,
            " ", input$habitat_type,
            " ", input$selected_wyt)

    active_streamgage_data <-
      get_data(streamgage_duration_rating_curves, package = "habistat") |>
      filter((station_id == input$streamgage_id) &
               (run == input$selected_run) &
               (habitat == input$habitat_type) &
               (wy_group == input$selected_wyt)) |>
      glimpse()

    if (nrow(active_streamgage_data) > 0) {

      (active_streamgage_data |> pull(data))[[1]] |>
        mutate(dhsi_selected = case_when(
          input$habitat_type == "spawning" ~ durhsi_spawning,
          active_comid_gradient_class == "Valley Lowland" ~ durhsi_rearing_vl,
          TRUE ~ durhsi_rearing_vf)) |>
        glimpse()

    } else {

      tibble(model_q = list(), dhsi_selected = list())

    }

  })

  output$streamgage_selector <- renderUI({

    active_comid_watershed_name <- (attr |>
      filter(comid == selected_point$comid) |>
      pull(watershed_level_3))[[1]]

    streamgage_options <- streamgages[[active_comid_watershed_name]]
    selectInput(inputId = "streamgage_id",
                label = "Select Gage for Duration Analysis",
                choices = setNames(names(streamgage_options), streamgage_options),
                selected = names(streamgage_options)[[1]])
  })

  duration_curve <- reactive({

    message("fsa")
    fsa <-
      predictions |>
      filter(comid == selected_point$comid) |>
      filter(habitat == input$habitat_type) |>
      transmute(flow_cfs, selected_wua = !!sym(input$wua_var)) |>
      glimpse()

    if (nrow(streamgage_drc()) > 0) {

      message("duration hsi curve")
      habistat::duration_apply_dhsi_to_fsa_curve(
        fsa = fsa,
        fsa_q = flow_cfs,
        fsa_wua = selected_wua,
        drc = streamgage_drc(),
        drc_q = model_q,
        drc_dhsi = dhsi_selected) |>
        glimpse()

    } else {

        tibble(q = list(), durwua = list())

    }
  })

}
