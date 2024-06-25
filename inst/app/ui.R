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
                     selectInput("active_flow", "Select Flow (cfs)", as.list(all_flows), selected=800),
                     selectInput("wua_var", "Select Variable", list("Scale-Dependent Model (prior baseflow removal)" = "wua_per_lf_pred_SD_a",
                                                                    "Scale-Dependent Model (post-model baseflow removal)" = "wua_per_lf_pred_SD_p",
                                                                    "Scale-Normalized Model (prior baseflow removal)" = "wua_per_lf_pred_SN_a",
                                                                    "Scale-Normalized Model (post-model baseflow removal)" = "wua_per_lf_pred_SN_p",
                                                                    "Actual" = "wua_per_lf_actual"))),
                 div(id = "comidDetails",
                     h2("Details for Selected Reach"),
                     div(id = "fsaPlot",
                         plotOutput("fsa_plot")),
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
