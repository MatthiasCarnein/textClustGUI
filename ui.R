library(shiny)
library(DT)
library(stopwords)
library(plotly)
library(shinyWidgets)


#plot the UI
ui <- shinyServer(fluidPage(
  sidebarLayout(
    sidebarPanel(
      titlePanel("Display", window="textClust"),
      tags$style(type = "text/css", ".irs-grid-text {font-size: 0px;}"), ## remove tick labels
      sliderTextInput("observation", "Observation", grid = TRUE, force_edges = TRUE, animate=animationOptions(interval = 2000), choices = c(0,0), selected = 0),
      prettyCheckbox("live", label = "Live", TRUE, icon = icon("check"), status="primary"),
      conditionalPanel("input.file_timeColumn!=null || input.streamType=='Server'",
                       prettyCheckbox("time_obs", label = "Show Timestamps", FALSE, icon = icon("check"), status="primary")
      ),
      radioGroupButtons(inputId="type", label="Type", choices=c("micro", "macro"), selected = "macro"),
      sliderTextInput("tokens","Tokens", grid = TRUE, force_edges = TRUE, choices=c(1:10,15,20, Inf), selected = "5"),
      
      tags$hr(),
      titlePanel("Stream"),
      radioGroupButtons(inputId="streamType", label="Source", choices=c("Server", "File", "Text"), selected = "Server", checkIcon=list(yes=icon("chevron-right"))),
      conditionalPanel("input.streamType=='Server'",
                       textInput("url", "URL", ""),
                       textInput("user", "User", ""),
                       passwordInput("password", "Password", ""),
                       numericInput("port", "Port", 27017),
                       actionButton("connect", "Connect", icon("database")),
                       actionButton("disconnect", "Disconnect", icon("ban")),
                       selectInput("streamSelect", "Select Stream", NULL),
                       sliderTextInput("refreshRate", "Refresh Rate (in sec)", grid = TRUE, force_edges = TRUE, choices = c(1:60, "None"), selected = "10")
      ),
      conditionalPanel("input.streamType=='Text'",
                       textAreaInput("streamInput", "Input Text", height = "100px")
      ),
      conditionalPanel("input.streamType=='File'",
                       fileInput("file", "Choose File", multiple = FALSE, accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv") 
                       )
      ),
      
      conditionalPanel("input.streamType=='File' || input.streamType=='Text'",
                       dropdown(
                         numericInput("file_textColumn", "Text Column", min = 1, val=1, step=1),
                         numericInput("file_timeColumn", "Time Column", min = 1, val=NULL, step=1),
                         textInput("file_separator", "Separator", value = "\t"),
                         conditionalPanel("input.file_timeColumn!=null",
                           textInput("file_timeFormat", "Time Format", value = "%Y-%m-%d %H:%M:%S"),
                           selectInput("file_timePrecision", "Time Precision", c("days", "hours", "minutes", "seconds"))),
                          prettyCheckbox("file_conversations", label = "Group By Conversations", FALSE, icon = icon("check"), status="primary"),
                         conditionalPanel("input.file_conversations",
                           numericInput("file_groupByColumn", "groupBy Column", min = 1, val=NA_integer_, step=1),
                           numericInput("file_parentTextColumn", "Parent Text Column", min = 1, val=NA_integer_, step=1),
                           numericInput("file_parentTimeColumn", "Parent Time Column", min = 1, val=NA_integer_, step=1)
                          ),
                         label = "Columns"
                       ),
                       sliderInput("stream_horizon", "Batch size", value = 100, min = 1, max=500, step=1),
                       actionButton("stream_start", "Start", icon("play")),
                       actionButton("stream_pause", "Pause", icon("pause")),
                       actionButton("stream_stop", "Stop", icon("stop")),
                       tags$hr(),
                       titlePanel("Algorithm"),
                       tags$div(title="Numer of clusters", numericInput("stream_setting_k", "k", value = 5, min = 1, step=1)),
                       tags$div(title="Radius threshold", numericInput("stream_setting_r", "r", min = 0, max = 2, value = .6, step=.1)),
                       tags$div(title="Fading factor", numericInput("stream_setting_lambda", "lambda", min = 0, max = 1, value = .002, step=.001)),
                       tags$div(title="Cleanup interval", numericInput("stream_setting_tgap", "tgap", min = 0, value = 100, step=10)),
                       tags$div(title="Minimum weight of micro clusters", numericInput("stream_setting_minWeight", "minWeight", min=0, value=3, step=.1)),
                       tags$div(title="Linkage method for reclustering", selectInput("stream_setting_linkage", "Linkage", c("single", "average", "complete"), selected = "complete")),
                       tags$div(title="Range of n-grams to calculate", sliderInput("stream_setting_ngrams", "N-Grams", min = 1, max = 5, value = c(1,2), step=1)),
                       tags$div(title="Fade tokens within clusters", prettyCheckbox("stream_setting_termFading", "Term Fading", TRUE, icon = icon("check"), status="primary")),
                       tags$div(title="Weight micro-clusters by their size during reclustering", prettyCheckbox("stream_setting_weightedReclustering", "Weighted Reclustering", TRUE, icon = icon("check"), status="primary")),
                       tags$div(title="Stopwords language and database",pickerInput(
                         inputId = "stream_setting_stopwords_sources", 
                         label = "Remove Stopwords", 
                         choices = lapply(sapply(stopwords_getsources(), function(x){
                           return(paste(x, stopwords_getlanguages(x), sep=" - "))
                         }), as.list), 
                         options = list(
                           `actions-box` = TRUE,
                           `deselect-all-text` = "None",
                           `select-all-text` = "All"
                         ),
                         multiple = TRUE,
                         selected = c("snowball - de", "snowball - en")
                       )),
                       fileInput("import", "Import", multiple = FALSE, accept=c(".RData"))
      ),
      conditionalPanel("input.streamType=='File' || input.streamType=='Text' || input.streamType=='Server'",
        downloadButton('export', "Export"))

    ),
    
    # main panel
    mainPanel(
      titlePanel("textClust"),
      # tab panel
      tabsetPanel(id="tabs",
                  tabPanel("N-Grams", 
                           dropdownButton(
                             prettyCheckbox("ngrams_showUnassigned", label = "Show Discarded Clusters", FALSE, icon = icon("check"), status="primary"),
                             downloadButton('downloadNgramData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           dataTableOutput("ngrams")),
                  tabPanel("Barplot", 
                           dropdownButton(
                             prettyCheckbox("barplot_showUnassigned", label = "Show Discarded Clusters", FALSE, icon = icon("check"), status="primary"),
                             sliderInput(inputId="barplot_numPlots", label="Number of Plots", value = c(1,10), min=1, max=10),
                             radioGroupButtons(inputId="barplot_plotType", label="", choices=c("Individual", "Grouped"), selected = "Individual"),
                             selectInput(inputId="barplot_order", label="Order", choices=c("Natural", "Weight", "Alphabetical"), selected = "Weight"),
                             downloadButton('downloadBarplot', "Download Plot"),
                             downloadButton('downloadBarplotData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           plotlyOutput("barplot",height = "800px")
                  ), 
                  tabPanel("MDS", 
                           dropdownButton(
                             prettyCheckbox("mds_showUnassigned", label = "Show Discarded Clusters", FALSE, icon = icon("check"), status="primary"),
                             radioGroupButtons(inputId="mds_plotType", label="", choices=c("Bubble", "Text"), selected = "Text"),
                             radioGroupButtons(inputId="mds_dimensionality", label="Dimensions", choices=c("2", "3"), selected = "2"),
                             downloadButton('downloadMDS', "Download Plot"),
                             downloadButton('downloadMDSData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           plotlyOutput("MDS",height = "800px")
                  ), 
                  tabPanel("t-SNE",
                           dropdownButton(
                             prettyCheckbox("tsne_showUnassigned", label = "Show Discarded Clusters", FALSE, icon = icon("check"), status="primary"),
                             radioGroupButtons(inputId="tsne_plotType", label="", choices=c("Bubble", "Text"), selected = "Text"),
                             radioGroupButtons(inputId="tsne_dimensionality", label="Dimensions", choices=c("2", "3"), selected = "2"),
                             sliderInput("tsne_iterations", "Max. Iterations", min = 1, max = 10000, value = 1000, step = 100),
                             sliderInput("tsne_perplexity", "Perplexity", min = 1, max = 50, value = 1, step = 1),
                             sliderInput("tsne_theta", "theta", min = 0, max = 1, value = .5, step = .05),
                             prettyCheckbox("tsne_pca", label = "PCA", FALSE, icon = icon("check"), status="primary"),
                             downloadButton('downloadTsne', "Download Plot"),
                             downloadButton('downloadTsneData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           plotlyOutput("tsne",height = "800px")
                  ),
                  tabPanel("Assignment", 
                           conditionalPanel("input.streamType=='File' || input.streamType=='Text' ||  input.streamType=='Server'",
                                            numericInput("selectCluster", "Cluster", value=1, min=1, max=1), 
                                            dataTableOutput("assignment")
                           )
                  ),
                  tabPanel("Sankey", 
                           dropdownButton(
                             prettyCheckbox("sankey_showUnassigned", label = "Show Discarded Clusters", FALSE, icon = icon("check"), status="primary"),
                             prettyCheckbox("sankey_microAndMacro", label = "Show Micro And Macro", FALSE, icon = icon("check"), status="primary"),
                             numericInput("sankey_trimTexts", "Trim Text Size", min=0, value = NA),
                             downloadButton('downloadSankey', "Download Plot"),
                             downloadButton('downloadSankeyData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           plotlyOutput("sankey")
                  ),
                  tabPanel("Dendrogram",
                           dropdownButton(
                             downloadButton('downloadDendrogram', "Download Plot"),
                             downloadButton('downloadDendrogramData', "Download Data"),
                             circle = TRUE,
                             icon = icon("gear"), 
                             size="sm",
                             tooltip = tooltipOptions(title = "Settings")
                           ),
                           plotOutput("dendrogram", height="1000px")
                  )
      )
    )
  )))

