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
                 div(id = "mapControls",
                     h2("Map Controls"),
                     radioButtons("flowline_scope", "Select Calculation Scope", choices=c("comid", "mainstem"), selected="comid", inline=T),
                     #selectInput("habitat_type", "Select Habitat Type", c("rearing", "spawning"), selected="rearing"),
                     radioButtons("habitat_type", "Select Habitat Type", choices=c("rearing", "spawning"), selected="rearing", inline=T),
                     #selectInput("active_flow", "Select Flow (cfs)", as.list(all_flows), selected=1000),
                     sliderInput("active_flow", "Select Flow (cfs) to show on map", min=min(all_flows), max=max(all_flows), step=100, value=1000),
                     #numericInput("active_flow", "Select Flow (cfs) to show on map", min=min(all_flows), max=max(all_flows), step=100, value=1000),
                     selectInput("wua_var", "Select Variable to show on map", list("Scale-Dependent Model" = "wua_per_lf_pred_SD",
                                                                    "Scale-Dependent Model (post-model baseflow removal)" = "wua_per_lf_pred_SD_ph_bfc_rm",
                                                                    "Scale-Normalized Model" = "wua_per_lf_pred_SN",
                                                                    "Scale-Normalized Model (post-model baseflow removal)" = "wua_per_lf_pred_SN_ph_bfc_rm",
                                                                    "Actual" = "wua_per_lf_actual"))),
                 div(id = "comidDetails",
                     h2("Details for Selected Reach"),
                     div(id = "fsaPlot",
                         plotOutput("fsa_plot"),
                         uiOutput("streamgage_selector"),
                         selectInput("selected_run", "Select Run", choices=c("fall", "late fall", "spring", "winter", "steelhead"), selected="fall"),
                         radioButtons("selected_wyt", "Select Water Year Type", choices=c("Dry", "Wet"), selected="Dry", inline=T)),
                     div(id = "predTable",
                         h3("Habitat Data"),
                         DT::DTOutput("pred_table")),
                     div(id = "attrTable",
                         h3("Flowline Attributes"),
                         DT::DTOutput("attr_table"))
                 ),
                 div(id = 'loading_radio', "Loading data, please wait...", style = "display: none;", class="sidebar-message")),
               mainPanel(
                 width = 6, # main width plus sidebar width should add to 12
                 shinyjs::useShinyjs(),  # Initialize shinyjs
                 shinycssloaders::withSpinner(leafletOutput("main_map"), hide.ui=F)
                 #shinycssloaders::withSpinner(leafgl::leafglOutput("main_map"))
                 )
               )
    )
  )
)
