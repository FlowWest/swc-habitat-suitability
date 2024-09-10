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
                     shinyWidgets::radioGroupButtons("flowline_scope", "Select Calculation Scope", choices=c("comid", "mainstem"), selected="comid"),
                     shinyWidgets::radioGroupButtons("habitat_type", "Select Habitat Type", choices=c("rearing", "spawning"), selected="rearing"),
                     shinyWidgets::radioGroupButtons("wua_var", "Select Calculation Method", list("Scale-Dependent" = "wua_per_lf_pred_SD",
                                                                              "Scale-Normalized" = "wua_per_lf_pred_SN",
                                                                              "Actual" = "wua_per_lf_actual")),
                     div(id = "mapControls",
                       div(id = "control_active_flow_slider",
                           shinyWidgets::sliderTextInput("active_flow", "Select Flow (cfs) to show on map", choices=all_flows_idx, selected=1000, hide_min_max=T),
                           style="display:inline-block; width:85%"),
                       div(id = "control_active_flow_apply",
                           actionButton("activeFlowApplyButton" ,"Apply"),
                           style="display:inline-block; width:10%; vertical-align: bottom;")
                     )),
                 div(id = "fsaPlot",
                     h3("Flow-to-Suitable-Area Plot"),
                     uiOutput("units_selector"),
                     uiOutput("clicked_item_heading"),
                     div(id = "fsaPlot",
                         shinycssloaders::withSpinner(plotOutput("fsa_plot"), hide.ui=F),
                         div(id = "durationOptions",
                            h4("Duration Analysis"),
                            uiOutput("streamgage_selector"),
                            selectInput("selected_run", "Select Run", choices=c("fall", "late fall", "spring", "winter", "steelhead"), selected="fall"),
                            shinyWidgets::radioGroupButtons("selected_wyt", "Select Water Year Type", choices=c("Dry", "Wet"), selected="Dry")),
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
