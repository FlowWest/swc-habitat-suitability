shinyUI(
  tagList(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    navbarPage(
      title = "HabiStat",
      id = "tabs",
      collapsible = TRUE,
      tabPanel("Interactive Map",
               sidebarPanel(
                 width = 6,
                 div(id = "mainControls",
                     radioButtons("flowline_scope", "Select Calculation Scope", choices=c("comid", "mainstem"), selected="comid", inline=T),
                     #selectInput("habitat_type", "Select Habitat Type", c("rearing", "spawning"), selected="rearing"),
                     radioButtons("habitat_type", "Select Habitat Type", choices=c("rearing", "spawning"), selected="rearing", inline=T),
                     ),
                 div(id = "mapControls",
                     h3("Map Controls"),
                     #selectInput("active_flow", "Select Flow (cfs)", as.list(all_flows), selected=1000),
                     #sliderInput("active_flow", "Select Flow (cfs) to show on map", min=min(all_flows), max=max(all_flows), step=100, value=1000),
                     shinyWidgets::sliderTextInput("active_flow", "Select Flow (cfs) to show on map", choices=as.character(all_flows_idx), selected="100", hide_min_max=T),
                     radioButtons("wua_var", "Select Calculation Method", list("Scale-Dependent" = "wua_per_lf_pred_SD",
                                                                              #"Scale-Dependent Model (post-model baseflow removal)" = "wua_per_lf_pred_SD_ph_bfc_rm",
                                                                              "Scale-Normalized" = "wua_per_lf_pred_SN",
                                                                              #"Scale-Normalized Model (post-model baseflow removal)" = "wua_per_lf_pred_SN_ph_bfc_rm",
                                                                              "Actual" = "wua_per_lf_actual"), inline=T)
                     ),
                 div(id = "fsaPlot",
                     h3("Flow-to-Suitable-Area Plot"),
                     uiOutput("clicked_item_heading"),
                     div(id = "fsaPlot",
                         shinycssloaders::withSpinner(plotOutput("fsa_plot"), hide.ui=F),
                         div(id = "durationOptions",
                            h4("Duration Analysis"),
                            uiOutput("streamgage_selector"),
                            selectInput("selected_run", "Select Run", choices=c("fall", "late fall", "spring", "winter", "steelhead"), selected="fall"),
                            radioButtons("selected_wyt", "Select Water Year Type", choices=c("Dry", "Wet"), selected="Dry", inline=T)),
                            shinycssloaders::withSpinner(plotOutput("dur_plot"), hide.ui=F),
                         )
                     ),
                 div(id = "detailTables",
                     div(id = "predTable",
                         h3("Habitat Data"),
                         #actionButton(inputId = "showPredTable",label = "Show/Hide"),
                         DT::DTOutput("pred_table")),
                     div(id = "attrTable",
                         h3("Flowline Attributes"),
                         #actionButton(inputId = "showAttrTable",label = "Show/Hide"),
                         DT::DTOutput("attr_table")
                         )
                     )
                 ),
               mainPanel(
                 width = 6, # main width plus sidebar width should add to 12
                 shinyjs::useShinyjs(),  # Initialize shinyjs
                 shinycssloaders::withSpinner(leafletOutput("main_map"), hide.ui=F),
                 #numericInput("active_flow", "Select Flow (cfs) to show on map", min=min(all_flows), max=max(all_flows), step=100, value=1000),
                 #shinycssloaders::withSpinner(leafgl::leafglOutput("main_map"))
                 )
               )
    )
  )
)
