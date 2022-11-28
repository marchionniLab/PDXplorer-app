
pcaInput <- function(id) {
  ns <- NS(id)

  plot_type <- radioButtons(ns("plot_type"), "Plot Type:", choices = c("ggplot"), inline = TRUE)

  title <- textInput(ns("title"), label = "Plot title", value = "")

  pc_x <- selectInput(ns("pc_x"), label = "PC shown on x", selected = "PC1", choices = paste0("PC", 1:20))

  pc_y <- selectInput(ns("pc_y"), label = "PC shown on y", selected = "PC2", choices = paste0("PC", 1:20))

  nr_genes <- numericInput(ns("nr_genes"), label = "Nr of (most variable) genes to use:", value = 500, min = 50, max = 20000)

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

  remove_batch_effect <- selectInput(ns("remove_batch_effect"), label = "Remove batch effect?", choices = NULL, selected = NULL, multiple = FALSE)

  remove_samples <- selectInput(ns("remove_samples"), label = "Remove sample by name:", choices = NULL, selected = NULL, multiple = TRUE)

  include_samples <- selectInput(ns("include_samples"), label = "Include sample by name:", choices = NULL, selected = NULL, multiple = TRUE)


  tagList(
    title,
    pc_x,
    pc_y,
    nr_genes,
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



pcaOutput <- function(id) {
  ns <- NS(id)

  tagList(

    ## pca
    shinydashboard::box(
      width = 12, title = "Principal component analysis",
      fluidRow(
        column(
          width = 3,
          checkboxInput(ns("sample_labels"),
            "Display sample labels?",
            value = FALSE
          )
        ),
        column(
          width = 3,
          checkboxInput(ns("remove_rug"),
            "Remove rug lines?",
            value = FALSE
          )
        )
      ),
      fluidRow(
        column(
          width = 3,
          selectInput(ns("cat_color_by"),
            label = "Color points by: ",
            choices = list("gene expression" = "geneSymbol", "metadata" = "metadata"),
            selected = "metadata", multiple = FALSE
          )
        ),
        column(
          width = 3,
          selectInput(ns("color_by"),
            label = "", choices = NULL,
            selected = NULL, multiple = FALSE
          )
        ),
        column(
          width = 3,
          selectInput(ns("cat_size_by"),
            label = "Point size by: ",
            choices = list("gene expression" = "geneSymbol", "metadata" = "metadata"),
            selected = "geneSymbol", multiple = FALSE
          )
        ),
        column(
          width = 3,
          selectInput(ns("size_by"),
            label = "",
            choices = NULL, selected = NULL,
            multiple = FALSE
          )
        )
      ),
      fluidRow(
        column(
          width = 3,
          selectInput(ns("shape_by"),
            label = "Point shape by:",
            choices = NULL, selected = NULL,
            multiple = FALSE
          )
        ),
        column(
          width = 3,
          selectInput(ns("encircle"),
            label = "Encircle samples by",
            selected = NULL, choices = NULL
          )
        ),
        column(
          width = 3,
          selectInput(ns("color_palette"),
            label = "Select color palette:",
            choices = c(names(ggsci_db)),
            selected = "auto", multiple = FALSE
          )
        ),
        column(
          width = 3,
          selectInput(ns("legend_position"),
            label = "Legend position:", choices = c("bottom", "top", "right", "left", "none" = "none"),
            selected = "bottom", multiple = FALSE
          )
        )
      ),
      plotOutput(ns("pca"), height = "800px"),
      #            plotOutput(ns("pca"), height="800px", click = clickOpts(id =ns("plot_click"))),
      #         verbatimTextOutput(ns("click_info")),
      ggDownloadUI(ns("ggpca")),
      collapsible = TRUE
    )
  )
}


## SERVER ######################################################################
pcaMod <- function(input, output, session, dds, metadata) {
  updateSelectInput(session, "shape_by", choices = c("none", names(metadata)), selected = "none")
  updateSelectInput(session, "encircle", choices = c("none", names(metadata)), selected = "none")
  updateSelectInput(session, "remove_samples", choices = c("", colnames(dds)), selected = "")
  updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected = "type")
  updateSelectInput(session, "include_samples", choices = c("", row.names(metadata)), selected = NULL)
  updateSelectInput(
    session,
    "remove_batch_effect",
    choices = c("none", "source", "batch", "sample_origin", "purity"),
    selected = NULL
  )


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

  observeEvent(input$cat_color_by, {
    if (input$cat_color_by == "geneSymbol") {
      updateSelectizeInput(session, "color_by", choices = c("none", row.names(dds)), selected = "none", server = TRUE)
    }
    if (input$cat_color_by == "metadata") {
      updateSelectizeInput(session, "color_by", choices = c("none", names(metadata)), selected = "none", server = TRUE)
    }
  })
  observeEvent(input$cat_size_by, {
    if (input$cat_size_by == "geneSymbol") {
      updateSelectizeInput(session, "size_by", choices = c("none", row.names(dds)), selected = "none", server = TRUE)
    }
    if (input$cat_size_by == "metadata") {
      updateSelectizeInput(session, "size_by", choices = c("none", names(metadata)), selected = "none", server = TRUE)
    }
  })



  pca_data <- reactive({
    shiny::validate(
      need(length(timer()$retain_samples) > 0 | length(timer()$include_samples) > 1,
        message = "Select at least 2 samples to get started."
      )
    )
    withProgress(message = "Running principal component analysis", {
      incProgress(0, detail = "This may take a while...")
      object <- dds[, timer()$final_list_samples_to_include]

      vst_data <- assay(object, "vst")

      if (input$remove_batch_effect != "none") {
        tryCatch(vst_data <- limma::removeBatchEffect(vst_data, batch = colData(object)[[input$remove_batch_effect]]), error = function(e) vst_data <- assay(object, "vst"))
      }

      pca_mat <- prep_PCA(vst_data, n_VarGenes = input$nr_genes, scale_features = TRUE)
      pcres <- prcomp(pca_mat)
      pccoords <- return_PCA_coords(pcres, cell_names = colnames(object), n_PCs = 20)
      return(pccoords$cells)
    })
  }) %>%
    bindCache(timer()$final_list_samples_to_include, input$nr_genes, input$remove_batch_effect)


  makePcaPlot <- function() {
    shiny::validate(
      need(!is.null(pca_data()),
        message = "Select samples to get started."
      )
    )

    shiny::validate(
      need(length(timer()$final_list_samples_to_include) > 0,
        message = "Select samples to get started."
      )
    )


    df_to_plot <- pca_data()

    if (input$color_by != "none") {
      if (length(unique(colData(dds[, timer()$final_list_samples_to_include])[[input$color_by]])) > length(ggsci_db[[input$color_palette]]$default)) {
        color_set <- colorRampPalette(ggsci_db[[input$color_palette]]$default)(length(unique(colData(dds[, timer()$final_list_samples_to_include])[[input$color_by]])))
      } else {
        color_set <- as.vector(ggsci_db[[input$color_palette]]$default[1:length(unique(colData(dds[, timer()$final_list_samples_to_include])[[input$color_by]]))])
      }
    }

    if (input$color_palette == "auto") {
      set_colors <- T
    } else {
      set_colors <- F
    }

    p <- plot_pca_deseq(dds[, timer()$final_list_samples_to_include],
      pca_results = df_to_plot,
      pc_x = input$pc_x, pc_y = input$pc_y,
      color_by = n2n(input$color_by),
      size_by = n2n(input$size_by),
      shape_by = n2n(input$shape_by),
      circle_by = n2n(input$encircle),
      remove_rug = input$remove_rug,
      set_colors = set_colors
    )

    # add labels if wanted
    if (input$sample_labels == TRUE) {
      p <- p + ggrepel::geom_label_repel(aes(label = sample_name, fill = NULL))
    }

    # add title
    p <- p + ggtitle(input$title) + theme(legend.position = input$legend_position)


    if (input$color_by != "none" & input$color_palette != "ggplot" & input$color_palette != "auto" & input$cat_color_by == "metadata") {
      p <- p + scale_color_manual(values = color_set)
    }

    p

    #  print(p)

    #    code <- callModule(ggDownload, "ggpca", reactive({p}))
  }


  callModule(ggDownload, "ggpca", reactive({
    makePcaPlot()
  }))



  #
  #   output$pca <- renderCachedPlot({
  #
  #     shiny::validate(
  #       need( !is.null(pca_data()),
  #             message = "Select samples to get started.")
  #     )
  #
  #     shiny::validate(
  #       need(length(timer()$final_list_samples_to_include) > 0,
  #            message = "Select samples to get started.")
  #     )
  #
  #
  #     df_to_plot <-  pca_data()
  #
  #     if(input$color_by != "none"){
  #       if ( length(unique(colData(dds[,timer()$final_list_samples_to_include])[[input$color_by]])) > length(ggsci_db[[input$color_palette]]$default) ){
  #         color_set = colorRampPalette(ggsci_db[[input$color_palette]]$default)(  length(unique(colData(dds[,timer()$final_list_samples_to_include])[[input$color_by]]))  )
  #       } else {
  #         color_set = as.vector( ggsci_db[[input$color_palette]]$default[1: length(unique(colData(dds[,timer()$final_list_samples_to_include])[[input$color_by]]))])
  #       }
  #     }
  #
  #     if (input$color_palette == "auto" ){ set_colors = T } else { set_colors = F }
  #
  #     p <- plot_pca_deseq(dds[,timer()$final_list_samples_to_include], pca_results =  df_to_plot,
  #                         pc_x = input$pc_x, pc_y = input$pc_y,
  #                         color_by = n2n(input$color_by),
  #                         size_by = n2n(input$size_by),
  #                         shape_by = n2n(input$shape_by),
  #                         circle_by = n2n(input$encircle),
  #                         remove_rug = input$remove_rug,
  #                         set_colors=set_colors)
  #
  #     # add labels if wanted
  #     if(input$sample_labels == TRUE){
  #       p <- p + ggrepel::geom_label_repel(aes(label=sample_name, fill = NULL))
  #     }
  #
  #     # add title
  #     p <- p + ggtitle(input$title) +  theme(legend.position=input$legend_position)
  #
  #
  #     if(input$color_by != "none" & input$color_palette != "ggplot" & input$color_palette != "auto"  & input$cat_color_by == "metadata"){
  #       p <- p + scale_color_manual(values=color_set)
  #     }
  #
  #     print(p)
  #
  #
  #   }, cacheKeyExpr = { list(pca_data(),
  #                            input$pc_x,
  #                            input$pc_y,
  #                            input$color_palette,
  #                            input$cat_color_by,
  #                            input$title,
  #                            input$legend_position,
  #                            input$remove_rug,
  #                            input$encircle,
  #                            input$shape_by,
  #                            input$size_by,
  #                            input$color_by,
  #                            timer()$final_list_samples_to_include,
  #                            input$sample_labels) })





  output$pca <- renderCachedPlot(
    {
      shiny::validate(
        need(!is.null(pca_data()),
          message = "Select samples to get started."
        )
      )

      shiny::validate(
        need(length(timer()$final_list_samples_to_include) > 0,
          message = "Select samples to get started."
        )
      )

      makePcaPlot()
    },
    cacheKeyExpr = {
      list(
        pca_data(),
        input$pc_x,
        input$pc_y,
        input$color_palette,
        input$cat_color_by,
        input$title,
        input$legend_position,
        input$remove_rug,
        input$encircle,
        input$shape_by,
        input$size_by,
        input$color_by,
        timer()$final_list_samples_to_include,
        input$sample_labels
      )
    }
  )
}
