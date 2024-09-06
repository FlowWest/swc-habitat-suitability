# TODO Duration gages should only be an option on mainstems
# Data not updating?
# (1) Check that the latest source hydraulic Rds outputs are present
# (2) Run model_cleaned.Rmd to re-export wua_hydraulic and wua_hydraulic_interp
# (3) Rebuild package and restart R session

function(input, output, session){

  # IDENTIFIERS FOR CLICKED ELEMENT --------------------------------------------

  most_recent_map_click <- reactiveValues(type = "none",
                                          lng = NULL,
                                          lat = NULL)

  clicked_item_label <- reactive({
    if (most_recent_map_click$type == "comid") {
      comid_name <- attr$gnis_name[[which(attr$comid==selected_point$comid)]]
      comid_watershed_name <- attr$watershed_level_3[[which(attr$comid==selected_point$comid)]]
      comid_name_suffix <- if_else(!is.na(comid_name),
                                   glue::glue(" ({comid_name})"),
                                   glue::glue(" (tributary of {comid_watershed_name})"))
      paste0("ComID ", selected_point$comid, comid_name_suffix)
    } else if (most_recent_map_click$type == "mainstem") {
      paste0(selected_mainstem$river_name, " Mainstem")
    } else if (most_recent_map_click$type == "watershed") {
      paste0(selected_watershed$watershed_name, " Watershed")
    }
  })

  output$clicked_item_heading <- renderUI({
    h4(clicked_item_label())
  })

  # REACTIVE FLOW FILTER (COMID) -----------------------------------------------

  message(paste(names(predictions), collapse = ", "))

  active_predictions <- reactive({
    # predictions_table |>
    #   filter(flow_cfs == input$active_flow) |>
    #   mutate(wua_per_lf = !!sym(input$wua_var))
    predictions |>
      filter(flow_idx == as.integer(input$active_flow)) |>
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
                                   #"Scale-Dependent Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SD_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model: ",
                                   round(wua_per_lf_pred_SN,2), " ft2/ft", "<br />",
                                   #"Scale-Normalized Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SN_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Actual: ", round(wua_per_lf_actual,2), " ft2/ft")) |>
      #filter(watershed_name == selected_watershed_name()) |>
      glimpse()
  })

  active_geom_mainstem <- reactive({
    mainstems |>
      inner_join(active_predictions_mainstem(),
                 by = join_by(river_cvpia)) |>
      mutate(wua_per_lf = !!sym(input$wua_var)) |>
      filter(!is.na(wua_per_lf)) |>
      mutate(object_label = paste0("<strong>", river_cvpia, "</strong>", "<br />",
                                   "Scale-Dependent Model: ",
                                   round(wua_per_lf_pred_SD,2), " ft2/ft", "<br />",
                                   #"Scale-Dependent Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SD_ph_bfc_rm,2), " ft2/ft", "<br />",
                                   "Scale-Normalized Model: ",
                                   round(wua_per_lf_pred_SN,2), " ft2/ft", "<br />")) |>
                                   #"Scale-Normalized Model (post-model baseflow removal): ",
                                   #round(wua_per_lf_pred_SN_ph_bfc_rm,2), " ft2/ft")
      #filter(watershed_name == selected_watershed_name()) |>
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
        most_recent_map_click$type <- "comid"
        selected_point$object_id <- input$main_map_shape_click$id
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
        selected_point$comid <- as.numeric(str_replace(selected_point$object_id, "comid_", ""))
      }
    }
  })

  # PREDICTION TABLES, ATTRIBUTES, PLOTS ---------------------------------------

  active_pred_table <- reactive({
    #if (!is.null(selected_point$object_id)) {
    if (most_recent_map_click$type == "comid") {
      active_predictions() |>
        filter(comid == selected_point$comid) |>
        select(where(is.numeric)) |>
        pivot_longer(everything()) |>
        mutate(across(value, function(x) round(x, 2)))
    } else if (most_recent_map_click$type == "watershed") {
      active_predictions_watershed() |>
        filter(watershed_level_3 == selected_watershed$watershed_id) |>
        select(where(is.numeric)) |>
        pivot_longer(everything()) |>
        mutate(across(value, function(x) round(x, 2)))
    }
  })

  output$pred_table <- DT::renderDT({
    if (most_recent_map_click$type == "none") {
      DT::datatable(active_pred_table(),
                    options = list(paging = F, searching = F))
    }
  })

  active_attr_table <- reactive({ # for COMID only
   #if (!is.null(selected_point$object_id) & (length(selected_point$comid)>0)) {
     # if (!is.null(selected_point$object_id)) {
     if (most_recent_map_click$type == "comid") {
     result <- attr |>
       filter(comid == selected_point$comid) |>
       select(where(is.numeric)) |>
       pivot_longer(everything())
       #gc()
     return(result)
   }
  })

  output$attr_table <- DT::renderDT({
     if (most_recent_map_click$type == "comid") {
        DT::datatable(active_attr_table(),
                      options = list(paging = F, searching = F))
      }
    })

  output$fsa_plot <- renderPlot({
    palette_linetypes <- c("Unscaled" = "solid",
                           "Duration Scaled" = "dashed")
    palette_colors <- c("Scale-Dependent" = "#6388b4",
                        "Scale-Normalized" = "#8cc2ca",
                        "Actual" = "#ffae34")

    if (most_recent_map_click$type == "comid") { #& (length(selected_point$comid)>0)) {
    predictions |>
      filter(comid == selected_point$comid) |>
      filter(habitat == input$habitat_type) |>
      ggplot(aes(x = flow_cfs)) +
      geom_line(aes(y = wua_per_lf_actual, color="Actual", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
      geom_line(aes(y = wua_per_lf_pred_SD, color="Scale-Dependent", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
      #geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
      geom_line(aes(y = wua_per_lf_pred_SN, color="Scale-Normalized", linetype = "Unscaled")) +#, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
     # geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
        input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
        input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
        input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
      #geom_hline(aes(yintercept = chan_width_ft)) + #, linetype="Channel Width (ft)")) +
      #geom_vline(aes(xintercept = baseflow_cfs)) +
      #geom_text(aes(x = 1, y = chan_width_ft, label = chan_width_ft)) +
      scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
      scale_y_continuous(limits = c(0, NA)) +
      #scale_y_continuous(trans = ihs, labels = scales::label_comma(), limits = c(0, NA)) +
      theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
      xlab("Flow (cfs)") + ylab("WUA (ft2) per linear ft") +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
      scale_linetype_manual(name = "Duration Analysis",
                            values = palette_linetypes)
    } else if (most_recent_map_click$type == "watershed") {
      #TODO: For watersheds, plot total acreage rather than WUA/LF
      predictions_watershed |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = wua_per_lf_pred_SD, color="Scale-Dependent", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = wua_per_lf_pred_SN, color="Scale-Normalized", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        #geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Analysis")) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab("WUA (ft2) per linear ft") +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
        scale_linetype_manual(name = "Duration Analysis",
                              values = palette_linetypes)
    } else if (most_recent_map_click$type == "mainstem") {
      predictions_mainstem |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = wua_per_lf_pred_SD, color="Scale-Dependent", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
       # geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = wua_per_lf_pred_SN, color="Scale-Normalized", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
          input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
          input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
          input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab("WUA (ft2) per linear ft") +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
        scale_linetype_manual(name = "Duration Analysis",
                              values = palette_linetypes)
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

  layer_flowlines <- function(m, show = TRUE, type = "comid") {
    if(type == "comid") {
    if(show) {
      # pal_limits <- c(min(active_geom()$wua_per_lf), max(active_geom()$wua_per_lf))
      message("plotting ", nrow(active_geom()), " flowlines")
      #m |> leafgl::addGlPolylines(data = active_geom(),layerId = ~object_id)
      m |> leaflet::addPolylines(data = active_geom(), # TODO switch to leafgl::addGlPolylines
                                 layerId = ~object_id,
                                 label = ~lapply(object_label, htmltools::HTML),
                                 color = ~pal(wua_per_lf, type=input$habitat_type), #pal(wua_per_lf),
                                 opacity = 1,
                                 weight = 2,
                                 options = leaflet::pathOptions(pane = "Flowlines"),
                                 group = "flowlines",
                                 highlightOptions = leaflet::highlightOptions(color = "#8B0000",
                                                                              weight = 3,
                                                                              bringToFront = TRUE)
      ) |>
        leaflet::addLegend(position = "bottomright",
                           colors = rev(flow_scale_colors[[input$habitat_type]]),
#                           colors = rev(pal(seq(pal_limits[[1]], pal_limits[[2]], (pal_limits[[2]]-pal_limits[[1]])/5))),
                           labels = lapply(paste("&ge;", rev(flow_scale_breaks[[input$habitat_type]])), htmltools::HTML),
#                           labels = c(round(pal_limits[[2]],2), rep("", 5-1), round(pal_limits[[1]],2)),
                           title = "Suitable Habitat Area (ft2) per linear ft",
                           layerId = "clegend")
      } else {
        m |> leaflet::removeShape(active_geom()$comid) |> leaflet::removeControl("clegend")
        # TODO switch to leaflet::removeGlPolylines
      }
    } else if (type == "mainstem") {
    if(show) {
      pal_limits <- c(min(active_geom_mainstem()$wua_per_lf), max(active_geom_mainstem()$wua_per_lf))
      message("plotting ", nrow(active_geom_mainstem()), " flowlines")
      m |> leaflet::addPolylines(data = active_geom_mainstem(),
                                 layerId = ~mainstem_id,
                                 label = ~lapply(object_label, htmltools::HTML),
                                 color = ~pal(wua_per_lf, type=input$habitat_type), #pal(wua_per_lf),
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
        m |> leaflet::removeShape(active_geom_mainstem()$object_id) |> leaflet::removeControl("clegend")
      }
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
      layer_flowlines(type = input$flowline_scope) #|>
    #leaflet::addLayersControl(baseGroups = c("watersheds", "none"),
    #                          overlayGroups = c("flowlines", "aoi"),
    #                          position = "bottomleft",
    #                          options = leaflet::layersControlOptions(collapsed = FALSE))
  })


# SELECTED WATERSHED OBSERVER ----------------------------------------------------

selected_watershed <- reactiveValues(object_id = NA,
                                     watershed_name = NA,
                                     lng = NA,
                                     lat = NA)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 10) == "watershed_") {
        most_recent_map_click$type <- "watershed"
        selected_watershed$object_id <- input$main_map_shape_click$id
        selected_watershed$watershed_name <-watersheds$watershed_level_3[[which(watersheds$watershed_id==input$main_map_shape_click$id)]]
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
      }
    }
  })

  selected_watershed_geom <- reactive({
    watersheds |>
      filter(watershed_id == selected_watershed$object_id) |>
      mutate(object_id = "active_watershed")
  })

  selected_watershed_bbox <- reactive({
    st_bbox(selected_watershed_geom())
  })

  selected_watershed_name <- reactive({
    if(!is.na(selected_watershed$object_id)) {
      watersheds$watershed_level_3[[which(watersheds$watershed_id==selected_watershed$object_id)]]
    } else {
      NA
    }
  })

  # SELECTED MAINSTEM OBSERVER ----------------------------------------------------

  selected_mainstem <- reactiveValues(object_id = NA,
                                       river_name = NA,
                                       lng = NA,
                                       lat = NA)

  observeEvent(input$main_map_shape_click, {
    cat(input$main_map_shape_click$id)
    if (!is.null(input$main_map_shape_click$id)) {
      if(substr(input$main_map_shape_click$id, 1, 9) == "mainstem_") {
        most_recent_map_click$type <- "mainstem"
        selected_mainstem$object_id <- input$main_map_shape_click$id
        selected_mainstem$river_name <-mainstems$mainstem_label[[which(mainstems$mainstem_id==input$main_map_shape_click$id)]]
        most_recent_map_click$lng <- input$main_map_shape_click$lng
        most_recent_map_click$lat <- input$main_map_shape_click$lat
      }
    }
  })

  active_mainstem_geom <- reactive({
    watersheds |>
      filter(mainstem_id == selected_mainstem$object_id) |>
      mutate(object_id = "active_mainstem")
  })

  # ACTIVE PLOT LOGIC ------------------------------------------------

  observe({

    proxy <- leaflet::leafletProxy("main_map") |>
      leaflet::removeShape("active_watershed")

    if (most_recent_map_click$type == "watershed") {
      # (nrow(selected_watershed_geom()) > 0) &
      proxy |>
        leaflet::addPolygons(data = selected_watershed_geom(),
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
        leaflet::flyToBounds(lng1 = selected_watershed_bbox()$xmin,
                             lat1 = selected_watershed_bbox()$ymin,
                             lng2 = selected_watershed_bbox()$xmax,
                             lat2 = selected_watershed_bbox()$ymax)
    } else {

      leaflet::leafletProxy("main_map") |>
        leaflet::removeShape("active_watershed")

    }

    if (most_recent_map_click$type == "comid") {
      proxy |>
        leaflet::addPolylines(data = active_geom(),
                             stroke = T,
                             weight = 2,
                             color = "red",
                             opacity = 1,
                             layerId = "active_comid",
                             label = ~lapply(object_label, htmltools::HTML))

    } else {

      leaflet::leafletProxy("main_map") |>
        leaflet::removeShape("active_comid")

    }

    if (most_recent_map_click$type == "mainstem") {
      proxy |>
        leaflet::addPolylines(data = active_geom_mainstem(),
                              stroke = T,
                              weight = 2,
                              color = "red",
                              opacity = 1,
                              layerId = "active_mainstem",
                              label = ~lapply(mainstem_label, htmltools::HTML))

    } else {

      leaflet::leafletProxy("main_map") |>
        leaflet::removeShape("active_mainstem")

    }


  })

  # REACTIVE FLOW FILTER (WATERSHED AND MAINSTEM) ------------------------------

  active_predictions_watershed <- reactive({
    predictions_watershed |>
      filter(flow_idx == as.integer(input$active_flow)) |>
      filter(habitat == input$habitat_type) |>
      glimpse()
  })

  active_predictions_mainstem <- reactive({
    predictions_mainstem |>
      filter(flow_idx == as.integer(input$active_flow)) |>
      filter(habitat == input$habitat_type) |>
      glimpse()
  })

  # DURATION ANALYSIS ----------------------------------------------------------

  streamgages <- get_data(streamgage_attr, package = "habistat") |>
    transmute(watershed, station_id,
              station_label = glue::glue("{str_to_upper(station_id)}: {str_to_upper(name)} ({min_wy}-{max_wy})")) |>
    nest(.by = watershed) |>
    mutate(data = map(data, function(x) deframe(x) |> as.list())) |>
    deframe()

  streamgage_locs <- get_data(streamgage_geom, package = "habistat") |>
    st_transform("+proj=longlat +datum=NAD83") |>
    st_set_crs("+proj=longlat +datum=WGS84") # for display purposes only

  streamgage_drc <- reactive({
    if (most_recent_map_click$type == "comid") {
    active_reach_gradient_class <- (attr |>
                                      filter(comid == selected_point$comid) |>
                                      pull(hqt_gradient_class))[[1]]
    } else if (most_recent_map_click$type == "mainstem") {
      # TODO: FIX (this is a placeholder using VF for all)
      active_reach_gradient_class <- "Valley Foothill"
    }

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
          active_reach_gradient_class == "Valley Lowland" ~ durhsi_rearing_vl,
          TRUE ~ durhsi_rearing_vf)) |>
        glimpse()

    } else {

      tibble(model_q = list(), dhsi_selected = list(), avg_max_days_inundated = list())

    }

  })

  output$streamgage_selector <- renderUI({

    active_reach_is_mainstem <- FALSE
    active_reach_watershed_name <- NA

    if (most_recent_map_click$type == "comid") {
      active_reach_watershed_name <- (attr |>
                                        filter(comid == selected_point$comid) |>
                                        pull(watershed_level_3))[[1]]
      active_reach_is_mainstem <- !is.na((attr |>
                                     filter(comid == selected_point$comid) |>
                                     pull(river_cvpia))[[1]])
    } else if (most_recent_map_click$type == "watershed"){
      active_reach_watershed_name <- selected_watershed$watershed_name
      active_reach_is_mainstem <- FALSE
    } else if (most_recent_map_click$type == "mainstem"){
      active_reach_watershed_name <- selected_mainstem$river_name
      active_reach_is_mainstem <- TRUE
    }

    if(active_reach_is_mainstem) {
      streamgage_options <- streamgages[[active_reach_watershed_name]]
      streamgage_options_geom <- streamgage_geom |> filter(station_id %in% names(streamgage_options))

      streamgage_options_selected <-
        streamgage_options_geom$station_id[[
          st_nearest_feature(st_point(c(most_recent_map_click$lng,
                                        most_recent_map_click$lat)),
                           streamgage_options_geom,
                           check_crs = FALSE)]]

      selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = setNames(names(streamgage_options), streamgage_options),
                  selected = streamgage_options_selected) # names(streamgage_options)[[1]])
    } else {
      selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = c())
    }
    # TODO Update this to default to the gage closest to the click point
  })

  duration_curve <- reactive({

    message("fsa")

    if (most_recent_map_click$type == "comid") {
      fsa <-
        predictions |>
        filter(comid == selected_point$comid) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(input$wua_var)) |>
        glimpse()
    } else if (most_recent_map_click$type == "mainstem") {
      fsa <-
        predictions_mainstem |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(input$wua_var)) |>
        glimpse()
    } else if (most_recent_map_click$type == "watershed") {
      fsa <-
        predictions_watershed |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(input$wua_var)) |>
        glimpse()
    } else {
      fsa <- tibble()
    }

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

  output$dur_plot <- renderPlot({
    if((nrow(streamgage_drc())>0) & (nrow(duration_curve()) > 0)) {
      bind_rows("days" = streamgage_drc() |> transmute(x = model_q, y = avg_max_days_inundated),
                "dhsi" = streamgage_drc() |> transmute(x = model_q, y = dhsi_selected),
                "wua" = bind_rows("wua" = duration_curve() |> transmute(x = q, y = wua),
                                  "durwua" = duration_curve() |> transmute(x = q, y = durwua),
                                  .id = "grp"),
                .id = "varname") |>
        bind_rows(tribble(~varname, ~grp, ~x, ~y,
                          "dhsi", "ref", -Inf, 1,
                          "dhsi", "ref", Inf,  1)) |>
        mutate(varname = varname |> factor(levels = c("days", "dhsi", "wua"),
                                           labels = c("Max Length of Period Exceeding Flow per WY-Season (days)",
                                                      "Duration Suitability Factor",
                                                      "Suitable Habitat Curve (ft2/ft)")),
               grp = grp |> coalesce("None") |>
                            factor(levels = c("None", "ref", "wua", "durwua"),
                                   labels = c("None", "ref", "Original", "Duration-Weighted"))) |>
        ggplot() +
        geom_step(aes(x = x, y = y, linetype = grp)) +
        facet_wrap(~varname, ncol = 1, scales="free_y") +
        scale_x_log10(breaks = scales::breaks_log(8), labels = scales::label_comma()) +
        scale_y_continuous(breaks = scales::breaks_extended(8), labels = scales::label_comma(), limits=c(0, NA)) +
        annotation_logticks(sides = "b") +
        ylab("") + xlab("Flow (cfs)") +
        theme(legend.position = "none",
              panel.grid.minor = element_blank()) +
        scale_linetype_manual(name = "",
                              values = c("None" = "solid",
                                         "ref" = "dashed",
                                         "Original" = "solid",
                                         "Duration-Weighted" = "dashed"))
    } else {
      ggplot()
    }

#    plt_days <-
#      ggplot() +
#      geom_step(data=streamgage_drc(), aes(x = model_q, y = avg_max_days_inundated)) +
#      xlab("Flow (cfs)") +
#      ylab("Days per WY") + scale_y_continuous(limits = c(0, 365)) +
#      ggtitle("Max Length of Period Exceeding Flow per WY-Season")
#
#    plt_dhsi <-
#      ggplot() +
#      geom_step(data=streamgage_drc(), aes(x = model_q, y = dhsi_selected)) +
#      xlab("Flow (cfs)") +
#      ylab("Weight") + scale_y_continuous(limits = c(0, 1)) +
#      ggtitle("Duration Suitability Factor")
#
#    plt_dwua <- ggplot() +
#      geom_line(data=duration_curve(), aes(x = q, y = wua, linetype = "Original")) +
#      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype = "Duration-Weighted")) +
#      xlab("Flow (cfs)") +
#      ylab("WUA (ft2/ft)") +
#      ggtitle("Suitable Habitat Area") +
#      scale_linetype_manual(name = "",
#                            values = c("Original" = "solid",
#                                       "Duration-Weighted" = "dashed"))
#
#    # plt_dwua <- set_dim(plt_dwua, get_dim(plt_dhsi))
#
#    (plt_days + plt_dhsi + plt_dwua) +
#      plot_layout(heights = c(2, 1, 2), axes = "collect", ncol = 1) &
#      theme(legend.position = "bottom") &
#      scale_x_log10(labels = scales::label_comma()) &
#      annotation_logticks(sides = "b")

  })

#   # INTERFACE
#
#   hide(id = "pred_table")
#   hide(id = "attr_table")
#
#   observeEvent(input$showPredTable,{
#     toggle(id = "pred_table",anim = T)
#   })
#
#   observeEvent(input$showAttrTable,{
#     toggle(id = "attr_table",anim = T)
#   })

}
