function(input, output, session){

  # REACTIVE FLOW FILTER -------------------------------------------------------

  active_predictions <- reactive({
    predictions_table |>
      filter(flow_cfs == input$active_flow) |>
      mutate(wua_per_lf = !!sym(input$wua_var))
  })

  active_geom <- reactive({
    flowlines_gcs |>
      inner_join(active_predictions(),
                 by=join_by(comid),
                 relationship="one-to-one") |>
      mutate(tooltip_html = paste0("<strong>", comid, " ", gnis_name, "</strong>", "<br />",
                                   "Scale-dependent model: ", round(wua_per_lf_pred_sd,2), " ft2/ft", "<br />",
                                   "Scale-normalized model: ", round(wua_per_lf_pred_si2,2), " ft2/ft", "<br />",
                                   "Two-step model: ", round(wua_per_lf_pred_sd2si,2), " ft2/ft",
                                   if_else(!is.na(wua_per_lf_actual),
                                           paste0("<br />", "Actual: ", round(wua_per_lf_actual,2), " ft2/ft"), "")))
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

  active_attr_table <- reactive({
    if (!is.null(selected_point$object_id)) {
      active_predictions() |>
        filter(comid == selected_point$comid) |>
        select(where(is.numeric), -wua_per_lf) |>
        pivot_longer(everything())
    }
  })

  output$attr_table <- DT::renderDT({
    if (!is.null(selected_point$object_id)) {
      DT::datatable(active_attr_table(),
                    options = list(paging = F, searching = F))
    }
  })

  output$fsa_plot <- renderPlot({
    if (!is.null(selected_point$object_id)) {
    predictions_table |>
      filter(comid == selected_point$comid) |>
      ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = wua_per_lf_pred_sd, color="Scale-Dependent Model")) +
        geom_line(aes(y = wua_per_lf_pred_si2, color="Scale-Normalized Model")) +
        geom_line(aes(y = wua_per_lf_pred_sd2si, color="Two-Step Model")) +
        geom_line(aes(y = wua_per_lf_actual, color="Actual")) +
        geom_hline(aes(yintercept = chan_width_ft, linetype="Channel Width (ft)")) +
        scale_x_log10(labels = scales::label_comma()) +
        #scale_y_continuous(trans = ihs, labels = scales::label_comma(), limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top") +
        xlab("Flow (cfs)") + ylab("WUA (ft2) per linear ft")
    } else {
      NULL
    }
  })

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
      pal_limits <- c(min(active_geom()$wua_per_lf), max(active_geom()$wua_per_lf))

      m |> leaflet::addPolylines(data = active_geom(),
                                 layerId = ~object_id,
                                 popup = ~tooltip_html,
                                 label = ~lapply(tooltip_html, htmltools::HTML),
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
    }
  }

  # LEAFLET RENDER -------------------------------------------------------------

  output$main_map <- renderLeaflet({
    shinyjs::showElement(id = 'loading_action')
    make_leaflet() |>
      layer_flowlines() #|>
      #leaflet::addLayersControl(baseGroups = c("watersheds", "none"),
      #                          overlayGroups = c("flowlines", "aoi"),
      #                          position = "bottomleft",
      #                          options = leaflet::layersControlOptions(collapsed = FALSE))
  })

}
