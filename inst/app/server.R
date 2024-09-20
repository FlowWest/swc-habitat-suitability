# TODO Duration gages should only be an option on mainstems
# Data not updating?
# (1) Check that the latest source hydraulic Rds outputs are present
# (2) Run model_cleaned.Rmd to re-export wua_hydraulic and wua_hydraulic_interp
# (3) Rebuild package and restart R session

function(input, output, session){

  # WUA VARIABLE UNITS ---------------------------------------------------------
  wua_var <- reactive({
    switch(input$wua_var,
           "wua_per_lf_pred_SD" = paste0(input$wua_units,"_pred_SD"),
           "wua_per_lf_pred_SN" = paste0(input$wua_units,"_pred_SN"),
           "wua_per_lf_actual" = paste0(input$wua_units,"_actual"))
    })

  wua_lab <- reactive({
    switch(input$wua_units,
           "wua_per_lf" = "Suitable Habitat Area (ft2) per linear ft",
           "wua_acres" = "Suitable Habitat Area (acres)")
  })

  wua_suf <- reactive({
    switch(input$wua_units,
           "wua_per_lf" = "ft2 / ft",
           "wua_acres" = "ac")
  })

  output$units_selector <- renderUI({
    shinyWidgets::radioGroupButtons("wua_units", "Select Habitat Units",
                 choices=c("ft2 per linear ft" = "wua_per_lf",
                           "total acres" = "wua_acres"),
                 selected= if (input$flowline_scope=="comid") "wua_per_lf" else "wua_acres",
                 disabled = is.null(input$main_map_shape_click$id))
  })

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
    #   filter(flow_cfs == active_map_params$flow) |>
    #   mutate(wua_per_lf = !!sym(input$wua_var))
    predictions |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
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
        select(any_of(names(var_names))) |>
        filter(comid == selected_point$comid) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    } else if (most_recent_map_click$type == "watershed") {
      active_predictions_watershed() |>
        select(any_of(names(var_names))) |>
        filter(watershed_level_3 == selected_watershed$watershed_id) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    } else if (most_recent_map_click$type == "mainstem") {
      active_predictions_mainstem() |>
        select(any_of(names(var_names))) |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
        pivot_longer(everything()) |>
        mutate(name = var_names[name])
    }
  })

  output$pred_table <- DT::renderDT({
    if (most_recent_map_click$type != "none") {
      DT::datatable(active_pred_table(),
                    options = list(paging = F, searching = F),
                    selection = "none")
    }
  })

  active_attr_table <- reactive({ # for COMID only
   #if (!is.null(selected_point$object_id) & (length(selected_point$comid)>0)) {
     # if (!is.null(selected_point$object_id)) {
     if (most_recent_map_click$type == "comid") {
     result <- attr |>
       filter(comid == selected_point$comid) |>
       mutate(across(where(is.numeric), function(x) signif(x, 3) |> as.character(x))) |>
       pivot_longer(everything()) |>
       mutate(name = var_names[name])
     #gc()
     return(result)
     }
  })

  output$attr_table <- DT::renderDT({
     if (most_recent_map_click$type == "comid") {
        DT::datatable(active_attr_table(),
                      options = list(paging = F, searching = F),
                      selection = "none")
      }
    })

  output$fsa_plot <- renderPlot({
    palette_linetypes <- c("Unscaled" = "solid",
                           "Duration Scaled" = "dashed")
    palette_colors <- c("Scale-Dependent" = "#6388b4",
                        "Scale-Normalized" = "#8cc2ca",
                        "Actual" = "#ffae34")
    #var_not_selected <- switch(input$wua_var,
    #                           "wua_per_lf_pred_SD"="wua_per_lf_pred_SN",
    #                           "wua_per_lf_pred_SN"="wua_per_lf_pred_SD")

    if (most_recent_map_click$type == "comid") { #& (length(selected_point$comid)>0)) {
    predictions |>
      filter(comid == selected_point$comid) |>
      filter(habitat == input$habitat_type) |>
      ggplot(aes(x = flow_cfs)) + #|> add_color_scale(type = input$habitat_type) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_actual")), color="Actual", linetype = "Unscaled")) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype = "Unscaled")) +
      geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype = "Unscaled")) +
      geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
        input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
        input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
        input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
      #geom_line(aes(y = !!sym(var_not_selected), linetype = "Unscaled"), color = "gray") +
      #geom_line(aes(y = !!sym(input$wua_var), linetype = "Unscaled"), linewidth=1) +
      #geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = durwua)) +
      #geom_hline(aes(yintercept = chan_width_ft)) + #, linetype="Channel Width (ft)")) +
      #geom_vline(aes(xintercept = baseflow_cfs)) +
      #geom_text(aes(x = 1, y = chan_width_ft, label = chan_width_ft)) +
      scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
      scale_y_continuous(limits = c(0, NA)) +
      #scale_y_continuous(trans = ihs, labels = scales::label_comma(), limits = c(0, NA)) +
      theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
      xlab("Flow (cfs)") + ylab(wua_lab()) +
      scale_color_manual(name = "Model Type",
                         values = palette_colors) +
      scale_linetype_manual(name = "Duration Analysis",
                            values = palette_linetypes)
    } else if (most_recent_map_click$type == "watershed") {
      predictions_watershed |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype = "Unscaled")) + #, linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        #geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Analysis")) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab(wua_lab()) +
        scale_color_manual(name = "Model Type",
                           values = palette_colors) +
        scale_linetype_manual(name = "Duration Analysis",
                              values = palette_linetypes)
    } else if (most_recent_map_click$type == "mainstem") {
      predictions_mainstem |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        filter(habitat == input$habitat_type) |>
        ggplot(aes(x = flow_cfs)) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SD")), color="Scale-Dependent", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
       # geom_line(aes(y = wua_per_lf_pred_SD_ph_bfc_rm, color="Scale-Dependent", linetype="Post-Model BFC Removal")) +
        geom_line(aes(y = !!sym(paste0(input$wua_units,"_pred_SN")), color="Scale-Normalized", linetype="Unscaled")) + # linetype=if_else(habitat=="rearing","Prior BFC Removal", "No BFC Removal"))) +
        #geom_line(aes(y = wua_per_lf_pred_SN_ph_bfc_rm, color="Scale-Normalized", linetype="Post-Model BFC Removal")) +
        geom_line(data=duration_curve(), aes(x = q, y = durwua, linetype="Duration Scaled", color = case_when(
          input$wua_var == "wua_per_lf_pred_SD" ~ "Scale-Dependent",
          input$wua_var == "wua_per_lf_pred_SN" ~ "Scale-Normalized",
          input$wua_var == "wua_per_lf_actual" ~ "Actual"))) +
        scale_x_log10(labels = scales::label_comma()) + annotation_logticks(sides = "b") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_minimal() + theme(panel.grid.minor = element_blank(), legend.position = "top", legend.box="vertical", text=element_text(size=21)) +
        xlab("Flow (cfs)") + ylab(wua_lab()) +
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
      leaflet::addMapPane("ValleyLowland", zIndex = 440) |>
      leaflet::addMapPane("Watersheds", zIndex = 445) |>
      leaflet::addMapPane("Flowlines", zIndex = 470) |>
      leaflet::addMapPane("Overlays", zIndex = 475) |>
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

    # first remove any existing flowlines
    m |> leaflet::removeShape(geom$object_id)
    m |> leaflet::removeShape(mainstems$object_id)

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
        m |> leaflet::removeShape(active_geom()$object_id) |> leaflet::removeControl("clegend")
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
                           colors = rev(flow_scale_colors[[input$habitat_type]]),
                           #                           colors = rev(pal(seq(pal_limits[[1]], pal_limits[[2]], (pal_limits[[2]]-pal_limits[[1]])/5))),
                           labels = lapply(paste("&ge;", rev(flow_scale_breaks[[input$habitat_type]])), htmltools::HTML),
                           #                           labels = c(round(pal_limits[[2]],2), rep("", 5-1), round(pal_limits[[1]],2)),
                           title = "Suitable Habitat Area (ft2) per linear ft",
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
                           weight = 1,
                           color = "red",
                           opacity = 0.5,
                           fill = T,
                           fillColor = "#FFFFFF",
                           fillOpacity = 0,
                           layerId = ~watershed_id,
                           group = "watersheds",
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
      leaflet::addPolygons(data = hqt, group = "HQT - Valley Lowland",
                           popup = "Habitat Quantification Tool Boundary",
                           color = "darkgrey",
                           weight = 0,
                           fillColor = "grey",
                           fillOpacity = 0.33,
                           options = leaflet::pathOptions(pane = "ValleyLowland")) |>
      addLayersControl(
        baseGroups = c("Basemap"),
        overlayGroups = c("flowlines", "watersheds", "HQT - Valley Lowland"), # TODO: unsure how to do watersheds right now..
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      hideGroup("HQT - Valley Lowland")
  })

  observe({
    proxy <- leaflet::leafletProxy("main_map", session = session)
    proxy |>
      layer_flowlines(type = input$flowline_scope)
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
        selected_mainstem$river_name <-mainstems$river_cvpia[[which(mainstems$mainstem_id==input$main_map_shape_click$id)]]
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

  # ACTIVE MAP LOGIC ------------------------------------------------

  observe({
    proxy <- leaflet::leafletProxy("main_map", session = session) |>
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
                             group = "watersheds", # this will make the active watershed also show/hide with the layer toggle
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
      clicked_polyline <- active_geom() |>
        filter(comid == selected_point$comid)

      bbox <- sf::st_bbox(clicked_polyline)

      proxy |>
        fitBounds(bbox$ymin, bbox$xmin, bbox$ymax, bbox$xmax) |> #TODO: get this to work...
        leaflet::addPolylines(data = clicked_polyline,
                             weight = 20,
                             color = "yellow",
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

  active_flow_applied <- eventReactive(input$activeFlowApplyButton, {
    input$active_flow
  })

  active_map_params <- reactiveValues(flow = 1000) # make sure this matches the value specified in the "selected" attribute on sliderTextInput in ui.R
  observe({
    active_map_params$flow <- active_flow_applied()
  })

  active_predictions_watershed <- reactive({
    predictions_watershed |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
      filter(habitat == input$habitat_type) |>
      glimpse()
  })

  active_predictions_mainstem <- reactive({
    predictions_mainstem |>
      filter(flow_idx == as.integer(active_map_params$flow)) |>
      filter(habitat == input$habitat_type) |>
      glimpse()
  })

  # DURATION ANALYSIS ----------------------------------------------------------

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

  # ACTIVE STREAMGAGE SELECTOR
  active_reach_info <- reactiveValues(is_mainstem = FALSE,
                                      watershed_name = NA,
                                      river_name = NA)

  observe({
    if (most_recent_map_click$type == "comid") {
      active_reach_info$watershed_name <- (attr |>
                                        filter(comid == selected_point$comid) |>
                                        pull(watershed_level_3))[[1]]
      active_reach_info$watershed_name <- (attr |>
                                             filter(comid == selected_point$comid) |>
                                             pull(river_cvpia))[[1]]
      active_reach_info$is_mainstem <- !is.na((attr |>
                                     filter(comid == selected_point$comid) |>
                                     pull(river_cvpia))[[1]])
    } else if (most_recent_map_click$type == "watershed"){
      active_reach_info$watershed_name <- selected_watershed$watershed_name
      active_reach_info$is_mainstem <- FALSE
    } else if (most_recent_map_click$type == "mainstem"){
      active_reach_info$river_name <- selected_mainstem$river_name
      active_reach_info$is_mainstem <- TRUE
    }
  })

  streamgage_options <- reactive({
    if(active_reach_info$is_mainstem) {
      streamgages[[active_reach_info$river_name]]
    } else {
      list()
    }
  })

  streamgage_options_geom <- reactive({
    message(paste(length(streamgage_options()), " streamgages found"))
    if(length(streamgage_options()) > 0){
      geom <- streamgage_pts |>
        filter(station_id %in% names(streamgage_options()))


      nearest_station_id <- geom$station_id[[
          st_nearest_feature(st_point(c(most_recent_map_click$lng,
                                        most_recent_map_click$lat)),
                             geom,
                             check_crs = FALSE)]]
      geom |> mutate(nearest = (station_id == nearest_station_id))
    } else {
      st_sf(st_sfc())
    }
  })

    output$streamgage_selector <- renderUI({
    if(active_reach_info$is_mainstem) {
      selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = setNames(names(streamgage_options()), streamgage_options()),
                  selected = streamgage_options_geom()$station_id[[which(streamgage_options_geom()$nearest)]])
    } else {
      selectInput(inputId = "streamgage_id",
                  label = "Select Gage for Duration Analysis",
                  choices = c())
    }
  })

    streamgage_options_geom_labelled <- reactive({
      if(length(streamgage_options()) > 0){
        streamgage_options_geom() |>
          mutate(selected = (station_id == input$streamgage_id))
      } else {
        st_sf(st_sfc())
      }
    })

  observe({
    streamgage_options_geom_labelled()
    proxy <- leaflet::leafletProxy("main_map")
    proxy |>
      leaflet::removeMarker(paste0("streamgage_", streamgage_pts$station_id))
    if (length(streamgage_options()) > 0) {
      proxy |>
        leaflet::addCircleMarkers(data = streamgage_options_geom_labelled(),
                                  layerId = ~paste0("streamgage_", station_id),
                                  group = "streamgages",
                                  label = ~station_label,
                                  color = ~if_else(selected, "#00A2E8", "#000000"),
                                  radius = 6,
                                  options = leaflet::markerOptions(pane = "Overlays"))
    }
  })

  observeEvent(input$main_map_marker_click, {
    cat(input$main_map_marker_click$id)
    if (!is.null(input$main_map_marker_click$id)) {
      if(substr(input$main_map_marker_click$id, 1, 11) == "streamgage_") {
        message("clicked streamgage")
        most_recent_map_click$lng <- input$main_map_marker_click$lng
        most_recent_map_click$lat <- input$main_map_marker_click$lat
      }
    }
  })

  # DURATION CURVE CALCULATION

  duration_curve <- reactive({

    message("fsa")

    if (most_recent_map_click$type == "comid") {
      fsa <-
        predictions |>
        filter(comid == selected_point$comid) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(wua_var())) |>
        glimpse()
    } else if (most_recent_map_click$type == "mainstem") {
      fsa <-
        predictions_mainstem |>
        filter(river_cvpia == selected_mainstem$river_name) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(wua_var())) |>
        glimpse()
    } else if (most_recent_map_click$type == "watershed") {
      fsa <-
        predictions_watershed |>
        filter(watershed_level_3 == selected_watershed$watershed_name) |>
        filter(habitat == input$habitat_type) |>
        transmute(flow_cfs, selected_wua = !!sym(wua_var())) |>
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
                                                      paste0("Suitable Habitat Curve (", wua_suf(),")"))),
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
