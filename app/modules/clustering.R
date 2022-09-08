
clusterInput <- function(id) {
  ns <- NS(id)


  title <- textInput(ns("title"), label = "Plot title", value = "")

  cor_method <- selectInput(
    inputId = ns("cor_method"),
    label = "Select correlation method",
    choices = c("pearson", "spearman"),
    selected = "spearman"
  )

  # https://github.com/UMMS-Biocore/debrowser_modular
  # https://onlinecourses.science.psu.edu/stat555/node/85
  dist_method <- selectInput(
    inputId = ns("dist_method"),
    label = "Select distance measure",
    choices = c("euclidean", "maximum", "manhattan", "canberra", "minkowski"),
    selected = "euclidean"
  )

  select_samples_by <- selectInput(
    inputId = ns("select_samples_by"),
    label = "Select samples by",
    choices = NULL,
    selected = NULL
  )

  selected_samples <- checkboxGroupInput(
    inputId = ns("selected_samples"),
    label = "",
    choices = NULL,
    selected = NULL
  )

  select_all <- actionLink(ns("select_all"), "Select all", style = "color:#18cc7b")
  deselect_all <- actionLink(ns("deselect_all"), "Deselect all", style = "color:#cc3c18")

  remove_samples <- selectInput(ns("remove_samples"), label = "Remove sample by name:", choices = NULL, selected = NULL, multiple = TRUE)

  include_samples <- selectInput(ns("include_samples"), label = "Include sample by name:", choices = NULL, selected = NULL, multiple = TRUE)

  remove_batch_effect <- selectInput(ns("remove_batch_effect"), label = "Remove batch effect?", choices = NULL, selected = NULL, multiple = FALSE)

  tagList(
    title,
    cor_method,
    dist_method,
    hr(),
    remove_batch_effect,
    hr(),
    h5("Include / exclude samples"),
    select_samples_by,
    tags$div(select_all, style = "display:inline-block; float:left"),
    tags$div(deselect_all, style = "display:inline-block; float:right"),
    br(),
    selected_samples,
    remove_samples,
    include_samples
  )
}



clusterOutput <- function(id) {
  ns <- NS(id)

  tagList(

    ## cluster
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Dendrogram",
        shinydashboard::box(
          width = 12, title = "Correlation dendrogram",
          fluidRow(
            column(
              width = 3,
              checkboxInput(ns("sample_labels"),
                "Display sample labels?",
                value = TRUE
              )
            ),
            column(
              width = 3,
              numericInput(ns("label_size"), label = "Label size: ", value = 2, min = 1, max = 10, step = 1)
            )
          ),
          fluidRow(
            column(
              width = 3,
              selectInput(ns("color_by"),
                label = "Color node by:", choices = NULL,
                selected = NULL, multiple = FALSE
              )
            ),
            column(
              width = 3,
              selectInput(ns("color_palette"),
                label = "Select color palette:",
                choices = c(names(ggsci_db)[-1]),
                selected = "ggplot", multiple = FALSE
              )
            ),
            column(
              width = 3,
              selectInput(ns("shape_by"),
                label = "Node shape by:",
                choices = NULL, selected = NULL,
                multiple = FALSE
              )
            ),
            column(
              width = 3,
              selectInput(ns("legend_position"),
                label = "Legend position:", choices = c("bottom", "top", "right", "left", "none" = "none"),
                selected = "bottom", multiple = FALSE
              )
            )
            # column(width = 3,
            #        selectInput(ns('color_label_by'), label = 'Color label by:', choices = NULL,
            #                    selected = NULL,multiple = FALSE))
          ),
          plotOutput(ns("dendro_dist"), height = "800px", click = clickOpts(id = ns("plot_click"))),
          ggDownloadUI(ns("ggdendro")),
          #     verbatimTextOutput(ns("click_info")),
          collapsible = FALSE
        )
      ),
      tabPanel(
        "Heatmap",
        shinydashboard::box(
          title = "Distance heatmap", width = 12,
          fluidRow(
            column(
              width = 3,
              checkboxInput(ns("sample_labels_hm"),
                "Display sample labels?",
                value = TRUE
              )
            ),
            column(
              width = 3,
              numericInput(ns("label_size_hm"), label = "Label size: ", value = 2, min = 1, max = 10, step = 1)
            )
          ),
          fluidRow(
            column(
              width = 3,
              selectInput(ns("column_annotation"),
                label = "Annotate row by:", choices = NULL,
                selected = NULL, multiple = TRUE
              )
            ),
            column(
              width = 3,
              selectInput(ns("color_palette_seq"),
                label = "Select color palette:",
                choices = row.names(subset(RColorBrewer:::brewer.pal.info, category == "div" | category == "seq")),
                selected = "RdYlBu", multiple = FALSE
              )
            ),
            column(
              width = 3,
              selectInput(ns("cutree_rows"),
                label = "Cut rows into N groups:", choices = seq(1, 20), selected = 1,
                multiple = FALSE
              )
            )
          ),
          plotOutput(ns("euc_dist_heatmap"), height = "800px"),
          ggDownloadUI(ns("ggheat")),
          collapsible = FALSE
        )
      )
    )
  )
}


## SERVER ######################################################################
clusterMod <- function(input, output, session, dds, metadata) {
  updateSelectInput(session, "shape_by", choices = c("none", names(metadata)), selected = "none")
  updateSelectInput(session, "color_by", choices = c("none", names(metadata)), selected = "none")
  updateSelectInput(session, "remove_samples", choices = c("", colnames(dds)), selected = "")
  updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected = "type")
  updateSelectInput(session, "column_annotation", choices = c("", names(metadata)), selected = "type")
  updateSelectInput(session, "include_samples", choices = c("", row.names(metadata)), selected = NULL)
  updateSelectInput(session, "remove_batch_effect", choices = c("none", "source", "batch", "patient", "purity"), selected = NULL)
  observeEvent(input$tabs, {
    if (input$tabs == "Dendrogram") {
      shinyjs::hide("dist_method")
      shinyjs::show("cor_method")
    } else if (input$tabs == "Heatmap") {
      shinyjs::show("dist_method")
      shinyjs::hide("cor_method")
    }
  })


  mydds <- reactiveValues(retain_samples = NULL, remove_samples = NULL, include_samples = NULL, final_list_samples_to_include = NULL)
  timer <- reactive(list(
    final_list_samples_to_include = mydds$final_list_samples_to_include,
    retain_samples = mydds$retain_samples,
    include_samples = mydds$include_samples
  )) %>% debounce(1000)

  observeEvent(input$select_samples_by,
    {
      updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected = unique(metadata[[input$select_samples_by]]), inline = TRUE)
    },
    ignoreInit = FALSE
  )
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected = unique(metadata[[input$select_samples_by]]), inline = TRUE)
  })

  observeEvent(input$deselect_all, {
    updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected = NULL, inline = TRUE)
  })
  observeEvent(input$selected_samples,
    {
      selected_column <- input$select_samples_by
      selected_groups <- input$selected_samples
      retain_samples <- row.names(metadata[metadata[[selected_column]] %in% selected_groups, ])

      columns_unselected <- setdiff(unique(metadata[[selected_column]]), selected_groups)
      rest_of_samples <- row.names(metadata[metadata[[selected_column]] %in% columns_unselected, ])
      if (length(input$selected_samples) > 0) {
        updateSelectInput(session, "remove_samples", choices = c(retain_samples), selected = input$remove_samples)
        updateSelectInput(session, "include_samples", choices = c(rest_of_samples), selected = input$include_samples)
        mydds$retain_samples <- retain_samples
      } else {
        mydds$retain_samples <- NULL
        updateSelectInput(session, "include_samples", choices = c("", row.names(metadata)), selected = input$include_samples)
      }
    },
    ignoreInit = FALSE,
    ignoreNULL = FALSE
  )
  observeEvent(input$remove_samples,
    {
      if (!is.null(input$remove_samples)) {
        outliersamples <- input$remove_samples
        mydds$remove_samples <- outliersamples
      } else {
        mydds$remove_samples <- NULL
      }
    },
    ignoreInit = FALSE,
    ignoreNULL = FALSE
  )
  observeEvent(input$include_samples,
    {
      if (!is.null(input$include_samples)) {
        includesamples <- input$include_samples
        mydds$include_samples <- includesamples
      } else {
        mydds$include_samples <- NULL
      }
    },
    ignoreInit = FALSE,
    ignoreNULL = FALSE
  )
  observe({
    mydds$final_list_samples_to_include <- setdiff(colnames(dds[, unique(c(mydds$retain_samples, mydds$include_samples))]), mydds$remove_samples)
  })



  dendro_data <- reactive({
    shiny::validate(
      need(length(timer()$retain_samples) > 0 | length(timer()$include_samples) > 0,
        message = "Select a sample to get started."
      )
    )

    withProgress(message = "Calculating correlation matrix for hierarchical cluster analysis", {
      incProgress(0)
      object <- dds[, timer()$final_list_samples_to_include]
      vst_data <- assay(object, "vst")

      if (input$remove_batch_effect != "none") {
        tryCatch(vst_data <- limma::removeBatchEffect(vst_data, batch = colData(object)[[input$remove_batch_effect]]), error = function(e) vst_data <- assay(object, "vst"))
      }

      dist <- as.dist(1 - cor(vst_data, method = input$cor_method))
      hc <- hclust(dist, method = "ward.D2")
      dendro <- as.dendrogram(hc)
      return(dendro)
    })
  }) %>%
    bindCache(timer()$final_list_samples_to_include, input$remove_batch_effect, input$cor_method)



  heatmapdata <- reactive({
    shiny::validate(
      need(length(timer()$retain_samples) > 0 | length(timer()$include_samples) > 0,
        message = "Select a sample to get started."
      )
    )

    withProgress(message = "Computing distance matrix for cluster analysis", {
      incProgress(0)
      object <- dds[, timer()$final_list_samples_to_include]
      vst_data <- assay(object, "vst")

      if (input$remove_batch_effect != "none") {
        tryCatch(vst_data <- limma::removeBatchEffect(vst_data, batch = colData(object)[[input$remove_batch_effect]]), error = function(e) vst_data <- assay(object, "vst"))
      }

      dist <- dist(t(vst_data), method = input$dist_method)
      dist_matrix <- as.matrix(dist)
      return(list(dist = dist, dist_matrix = dist_matrix))
    })
  }) %>%
    bindCache(timer()$final_list_samples_to_include, input$remove_batch_effect, input$dist_method)



  makeDendroPlot <- function() {
    shiny::validate(
      need(length(timer()$final_list_samples_to_include) > 1,
        message = "Select at least two samples to get started."
      )
    )


    if (input$color_by != "none") {
      colorby <- input$color_by
    } else {
      colorby <- NULL
    }
    if (input$shape_by != "none") {
      shapeby <- input$shape_by
    } else {
      shapeby <- NULL
    }

    hcd <- dendro_data()

    ddata_x <- ggdendro::dendro_data(hcd)
    P <- ggplot(ggdendro::segment(ddata_x)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend))
    labs <- ggdendro::label(ddata_x)
    ymax <- max(ddata_x$segments$yend)
    labelspace <- 0.2
    total_axis_size <- ymax * (1 / (1 - labelspace))
    labs <- cbind(labs, metadata[match(labs$label, row.names(metadata)), ])

    if (input$color_by != "none") {
      if (length(unique(labs[[colorby]])) > length(ggsci_db[[input$color_palette]]$default)) {
        color_set <- colorRampPalette(ggsci_db[[input$color_palette]]$default)(length(unique(labs[[colorby]])))
      } else {
        color_set <- as.vector(ggsci_db[[input$color_palette]]$default[1:length(unique(labs[[colorby]]))])
      }
    }



    if (input$sample_labels) {
      P <- P + geom_text(data = labs, angle = 90, hjust = 1, size = rel(input$label_size), aes_string(label = "label", x = "x", y = -(ymax / 40), colour = colorby), show.legend = F)
    }

    P <- P + ggdendro::theme_dendro() + ylim(-(total_axis_size * labelspace), ymax) #+ scale_color_manual(name = "patient", values = palette)

    if (input$shape_by != "none") {
      shapes <- rep(15:20, 10)[1:length(unique(labs[[shapeby]]))]
    } else {
      shapes <- 1
    }

    P <- P + geom_point(data = labs, aes_string(x = "x", y = 0, colour = colorby, shape = shapeby), size = 2.5) + scale_shape_manual(values = shapes, name = shapeby) + theme(title = element_text(size = rel(1)), legend.text = element_text(size = rel(1)))
    P <- P + ggtitle(input$title) + theme(plot.title = element_text(size = 16), legend.position = input$legend_position)

    if (input$color_palette != "ggplot" & input$color_by != "none") {
      P <- P + scale_color_manual(values = color_set)
    }

    P
  }

  output$dendro_dist <- renderCachedPlot(
    {
      shiny::validate(
        need(length(timer()$final_list_samples_to_include) > 1,
          message = "Select at least two samples to get started."
        )
      )

      makeDendroPlot()
    },
    cacheKeyExpr = {
      list(
        dendro_data(),
        input$title,
        input$color_by,
        input$shape_by,
        input$color_palette,
        input$sample_labels,
        timer()$final_list_samples_to_include,
        input$legend_position,
        input$remove_batch_effect,
        input$cor_method,
        input$label_size
      )
    }
  )



  makeHeatPlot <- function() {
    shiny::validate(
      need(length(timer()$final_list_samples_to_include) > 1,
        message = "Select at least two samples to get started."
      )
    )

    sampleDists <- heatmapdata()$dist
    sampleDistMatrix <- heatmapdata()$dist_matrix

    colnames(sampleDistMatrix) <- NULL

    if (length(input$column_annotation) > 0) {
      expgroups <- as.data.frame(colData(dds)[, input$column_annotation])
      rownames(expgroups) <- colnames(dds)
      colnames(expgroups) <- input$column_annotation
    } else {
      expgroups <- NULL
    }

    pheatmap::pheatmap(sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name = input$color_palette_seq)))(100),
      cutree_rows = input$cutree_rows,
      main = input$title,
      annotation_row = expgroups,
      fontsize_row = input$label_size_hm,
      show_rownames = input$sample_labels_hm
    )
  }

  output$euc_dist_heatmap <- renderCachedPlot(
    {
      shiny::validate(
        need(length(timer()$final_list_samples_to_include) > 1,
          message = "Select at least two samples to get started."
        )
      )

      makeHeatPlot()
    },
    cacheKeyExpr = {
      list(
        heatmapdata(),
        input$dist_method,
        input$column_annotation,
        input$remove_batch_effect,
        input$title,
        input$sample_labels_hm,
        input$label_size_hm,
        input$cutree_rows,
        input$color_palette_seq,
        timer()$final_list_samples_to_include
      )
    }
  )


  callModule(ggDownload, "ggheat", reactive({
    makeHeatPlot()
  }))
  callModule(ggDownload, "ggdendro", reactive({
    makeDendroPlot()
  }))
}
