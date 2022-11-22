library(ggplot2)
library(ABCutilities)
library(shiny)
library(shinythemes)
library(shinyjs)
library(data.table)
library(scales)
library(dplyr)
library(grid)
library(gridExtra)
library(codeModules)
library(ggdendro)

source("modules/pca.R")
source("modules/clustering.R")
source("modules/gene_explorer.R")
source("modules/fusions.R")
source("modules/selectSamples.R")
source("modules/ssgsea.R")
source("modules/differential.R")
source("modules/markers.R")

invisible(
  sapply(
    list.files(path = "R/", pattern = "\\.R$", full.names = TRUE),
    source
  )
)
shiny::shinyOptions(
  cache = cachem::cache_disk("./pdxrna-app-cache")
)


dds <- reactiveValues(
  dds = NULL,
  metadata = NULL,
  fusions = NULL,
  gsva = NULL,
  ssgsea = NULL
)

rev_num <- "0.1.0"

# Define UI for application that draws a histogram
ui <- function(request) {
  fluidPage(
    # theme = shinytheme("simplex"), useShinyjs(), # spacelab
    theme = bslib::bs_theme(
      version = 5,
      bootswatch = "simplex"
    ),
    tags$head(
      tags$style(
        HTML("
              shiny-output-error-validation {
              color: red;
              font-size: 14px;
              }
        ")
      )
    ),
    tags$style(type = "text/css", ".navbar-nav {padding-top: 8px;}"),
    # TODO: @luciorq Move CSS, JS, and HTML to relevant dirs
    # tags$head(
    #  HTML('<script defer src="https://use.fontawesome.com/releases/v5.0.12/js/all.js" integrity="sha384-Voup2lBiiyZYkRto2XWqbzxHXwzcm4A5RfdfG6466bu5LqjwwrjXCMBQBLMWh7qR" crossorigin="anonymous"></script>')
    # ),
    # added  https://github.com/rstudio/shinydashboard
    includeCSS(path = "AdminLTE.css"),
    #  https://github.com/rstudio/shinydashboard
    includeCSS(path = "shinydashboard.css"),
    #  https://github.com/rstudio/shinydashboard/blob/master/srcjs/AdminLTE/app.js
    includeScript(path = "app.js"),

    # Bookmark Button
    shiny::bookmarkButton(),
    navbarPage(
      title = HTML(
        paste0("PDXplorer RNA <br> <font size='0.5'>", rev_num, "</font>")
      ),
      windowTitle = "PDXplorer-RNA",
      id = "tabs",
      collapsible = TRUE,
      tabPanel("PCA",
        icon = icon("cube"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            pcaInput("pca")
          ),
          mainPanel(
            width = 9,
            pcaOutput("pca")
          )
        )
      ),
      tabPanel("Clustering",
        icon = icon("sitemap"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            clusterInput("cls")
          ),
          mainPanel(
            width = 9,
            clusterOutput("cls")
          )
        )
      ),
      tabPanel("Gene exploration",
        icon = icon("chart-bar"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            geInput("ge")
          ),
          mainPanel(
            width = 9,
            geOutput("ge")
          )
        )
      ),
      tabPanel("ssGSEA",
        icon = icon("arrow-up-wide-short"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            gseaInput("gsea")
          ),
          mainPanel(
            width = 9,
            gseaOutput("gsea")
          )
        )
      ),
      tabPanel("Differential genes",
        icon = icon("calculator"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            diffInput("diff")
          ),
          mainPanel(
            width = 9,
            diffOutput("diff")
          )
        )
      ),
      tabPanel("Marker genes",
        icon = icon("thumbtack"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            markersInput("markers")
          ),
          mainPanel(
            width = 9,
            markersOutput("markers")
          )
        )
      ),
      tabPanel("Fusion explorer",
        icon = icon("connectdevelop"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            fusionInput("fus")
          ),
          mainPanel(
            width = 9,
            fusionOutput("fus")
          )
        )
      )
    )
  )
}

server <- function(input, output, session) {
  RV <- reactiveValues(
    pca_ui_flg = FALSE,
    cls_ui_flg = FALSE,
    ge_ui_flg = FALSE,
    fus_ui_flg = FALSE,
    diff_ui_flg = FALSE,
    marker_ui_flg = FALSE,
    ssgsea_ui_flag = FALSE
  )

  observeEvent(input$tabs, {
    if (grepl("PCA", input$tabs)) {
      if (RV$pca_ui_flg) {
        return()
      }

      # TODO: @luciorq Use `ShinyCSSLoaders` while waiting for modules
      withProgress(message = "Loading dataset", {
        incProgress(0, detail = "This may take a while...")

        # TODO: @luciorq Move data (RDS) to a database or .parquet files
        dds$dds <- readRDS("../data/dds_all_nov2022.rds")
        colData(dds$dds)$sampleName <- make.names(
          as.character(colData(dds$dds)$sampleName)
        )
        colData(dds$dds)$type_source <- paste0(
          colData(dds$dds)$type, "_", colData(dds$dds)$source
        )
        # dds$dds$type_design <- make.names(dds$dds$type_design)
        dds$metadata <- data.frame(
          colData(dds$dds)
        )[, c("patient", "passage", "source", "type", "purity"), ]
        dds$metadata$type_source <- paste0(
          dds$metadata$type, "_", dds$metadata$source
        )

        # TODO: @luciorq Move from `callModule` to `moduleServer`
        callModule(
          pcaMod, "pca",
          dds = dds$dds,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # Gene Exploration
        callModule(
          geMod,
          "ge",
          dds = dds$dds,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # Clustering
        callModule(
          clusterMod,
          "cls",
          dds = dds$dds,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # Marker genes
        callModule(
          markersMod,
          "markers",
          dds = dds$dds,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # Fusion explorer
        dds$fusions <- readRDS("../data/fusions_2019.RDS")
        dds$fusions$summary$sample <- make.names(dds$fusions$summary$sample)
        dds$fusions$list$sample <- make.names(dds$fusions$list$sample)
        callModule(
          fusionMod,
          "fus",
          fusions = dds$fusions,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # ssGSEA
        dds$gsva <- readRDS("../data/gsva.RDS")
        dds$ssgsea <- readRDS("../data/ssgsea.RDS")
        callModule(
          gseaMod,
          "gsea",
          dds = dds$dds,
          gsva = dds$gsva,
          ssgsea = dds$ssgsea,
          metadata = dds$metadata
        )
        incProgress(1 / 7)

        # Differential genes
        dds$fit <- readRDS("../data/fit_all_march2021.RDS")
        dds$fit_pdx <- readRDS("../data/fit_pdxOnly_march2021.RDS")
        dds$fit_primary <- readRDS("../data/fit_primaryOnly_march2021.RDS")

        callModule(
          diffMod,
          "diff",
          dds = dds$dds,
          metadata = dds$metadata,
          fit = dds$fit,
          fit_pdx = dds$fit_pdx,
          fit_primary = dds$fit_primary
        )
        incProgress(1 / 7)
      })
      RV$pca_ui_flg <- TRUE
    }
  })
}

# Run the application
shiny::enableBookmarking("url")
shiny::shinyApp(
  ui = ui,
  server = server,
  enableBookmarking = "url"
)
