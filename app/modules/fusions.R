
fusionInput <- function(id) {

  ns <- NS(id)


  title <- textInput(ns("title"), label = "Plot title", value = "Number of fusions detected")

  facet_by_row <- selectInput(inputId=ns("facet_by_row"),
                              label="Facet plot row-wise by",
                              choices = NULL,
                              selected = NULL)

  facet_by_col <- selectInput(inputId=ns("facet_by_col"),
                              label="Facet plot column-wise by",
                              choices = NULL,
                              selected = NULL)

  select_samples_by <- selectInput(inputId=ns("select_samples_by"),
                                   label="Select samples by",
                                   choices = NULL,
                                   selected = NULL)

  selected_samples <- checkboxGroupInput(inputId=ns("selected_samples"),
                                         label="",
                                         choices = NULL,
                                         selected = NULL)

  selected_fusion <- selectizeInput(inputId=ns("selected_fusion"),
                                         label="Display samples with selected fusion:",
                                         choices = NULL,
                                         selected = NULL,
                                        multiple=TRUE)

  select_all <- actionLink(ns("select_all"), "Select all boxes")

  remove_samples <-selectInput(ns('remove_samples'), label = "Remove sample by name:", choices = NULL, selected=NULL, multiple = TRUE)

  include_samples <- selectInput(ns('include_samples'), label = "Include sample by name:", choices = NULL, selected=NULL, multiple = TRUE)

  set_size <- numericInput(ns('set_size'), label = "Show sets with at least n members:", min = 1, max = 100, value = 5)


   tagList(title,
           facet_by_row,
           facet_by_col,
          hr(),
          uiOutput(ns('ui_header')),
       #   h5("Include / exclude samples"),
          select_samples_by,
          select_all,
          selected_samples,
          remove_samples,
          include_samples,
         set_size,
          hr(),
          selected_fusion)

}



fusionOutput <- function(id) {
  ns <- NS(id)

  tagList(

    ## cluster
    tabsetPanel(id = ns("tabs"),
      tabPanel("Summary",
                 shinydashboard::box( width = 12, title = "Overview",
                                       fluidRow(column(width = 3,
                                                       checkboxInput(ns('sample_labels'),
                                                                     "Display sample labels?",
                                                                     value = TRUE)),
                                                column(width = 3,
                                                       checkboxInput(ns('percentage'),
                                                                     "Display as percentage?",
                                                                     value = FALSE)) ),
                                       fluidRow(
                                      #   column(width = 3,
                                      #          selectInput(ns('color_by'), label = 'Color node by:', choices = NULL,
                                      #                      selected = NULL,multiple = FALSE)),
                                        column(width = 3,
                                               selectInput(ns('color_palette'),
                                                           label = 'Select color palette:',
                                                           choices = c(names(ggsci_db)[-1]),
                                                           selected = "npg" ,multiple = FALSE)),

                                               column(width = 3,
                                                      selectInput(ns('legend_position'), label = 'Legend position:', choices = c("bottom", "top", "right", "left", "none" = "none"),
                                                                  selected = "bottom",multiple = FALSE))
                                               # column(width = 3,
                                               #        selectInput(ns('color_label_by'), label = 'Color label by:', choices = NULL,
                                               #                    selected = NULL,multiple = FALSE))
                                      ),
                                      plotOutput(ns("overview"), height="1000px", click = clickOpts(id =ns("plot_click"))),
                                      #     verbatimTextOutput(ns("click_info")),
                                      collapsible = FALSE)
      ),

      tabPanel("Tables",
                 shinydashboard::box(title = "Tables", width = 12,
                                     DT::dataTableOutput(ns('fusion_table'), height="1000px"),
                   collapsible = FALSE

               )
      ),

      tabPanel("Intersections",
               tabsetPanel(id = ns("inter_tabs"),
                           tabPanel("empty", "This is the hello tab")
               ),
               tags$div(id = ns('placeholder')),
               uiOutput(ns("inter_tables"))
      ),

      tabPanel("Circos plots",
                shinydashboard::box( width = 12, title = "Circos plot",
                                     fluidRow(
                                       column(width = 5,
                                              selectInput(ns('use_fusions_from'),
                                                          label = 'Select data from:',
                                                          choices = c("STAR Fusion" = "star.fusion_final", "Fusion Catcher" = "fusion.catcher_final"),
                                                          selected = "star.fusion_final",multiple = FALSE))
                                     ),
                                     plotOutput(ns("circos_plot"), height="800px", width="800px", click = clickOpts(id =ns("plot_click")))
                )

       )
  ))
}


## SERVER ######################################################################
fusionMod <- function(input, output, session, fusions, metadata) {


  updateSelectInput(session, "remove_samples", choices=c("",row.names(metadata)), selected="")
  updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected="type")
  updateSelectInput(session, "columun_annotation", choices = c("", names(metadata)), selected="")
  updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=NULL)
  updateSelectInput(session, "facet_by_row", choices = c(none = ".", names(metadata)), selected=".")
  updateSelectInput(session, "facet_by_col", choices = c(none = ".", names(metadata)), selected=".")

  data <- reactiveValues(data = NULL, retain_samples = NULL, remove_samples = NULL, include_samples = NULL, redraw = TRUE)

  hideTab(inputId = "inter_tabs", target = "empty")

  observeEvent(input$tabs, {
    if (input$tabs ==  "Summary"){
      output$ui_header <- renderUI({
        h5("Include / exclude samples")
      })
      shinyjs::enable("selected_samples")
      shinyjs::hide("selected_fusion")
      shinyjs::show("title")
      shinyjs::show("facet_by_col")
      shinyjs::show("facet_by_row")
      shinyjs::show("remove_samples")
      shinyjs::show("include_samples")
      shinyjs::show("selected_samples")
      shinyjs::show("select_all")
      shinyjs::hide("set_size")
      updateTextInput(session, "title",label = "Plot title", value ="Number of fusions detected")

    } else if (input$tabs ==  "Tables"){
      output$ui_header <- renderUI({
        h5("Include / exclude samples")
      })
      shinyjs::enable("selected_samples")
      shinyjs::show("selected_fusion")
      shinyjs::hide("title")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("facet_by_row")
      shinyjs::show("remove_samples")
      shinyjs::show("include_samples")
      shinyjs::show("selected_samples")
      shinyjs::show("select_all")
      shinyjs::hide("set_size")

    } else if (input$tabs == "Intersections"){
      output$ui_header <- renderUI({
        h5("Find intersections by: ")
      })
      shinyjs::hide("selected_fusion")
      shinyjs::hide("title")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("remove_samples")
      shinyjs::hide("include_samples")
      updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= unique(metadata[[input$select_samples_by]]), inline=TRUE )
      shinyjs::hide("selected_samples")
      shinyjs::hide("select_all")
      shinyjs::show("set_size")
    } else if (input$tabs ==  "Circos plots"){
      output$ui_header <- renderUI({
        h5("Include / exclude samples")
      })
      shinyjs::enable("selected_samples")
      shinyjs::hide("selected_fusion")
      shinyjs::show("title")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("facet_by_row")
      shinyjs::show("remove_samples")
      shinyjs::show("include_samples")
      shinyjs::show("selected_samples")
      shinyjs::show("select_all")
      shinyjs::hide("set_size")
      updateTextInput(session, "title",label = "Plot title", value = "Fusions for selected sample(s)")
    }


  })



  observeEvent(input$select_samples_by, {
    updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected=NULL, inline=TRUE )
  }, ignoreInit = FALSE)


  observeEvent(input$select_all,{
    updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= unique(metadata[[input$select_samples_by]]), inline=TRUE )
  })

  observeEvent(input$selected_samples, {
    selected_column <- input$select_samples_by
    selected_groups <- input$selected_samples
    retain_samples <- row.names(metadata[metadata[[selected_column]] %in% selected_groups,])
    columns_unselected <- setdiff(unique(metadata[[selected_column]]) , selected_groups)
    rest_of_samples <- row.names( metadata[metadata[[selected_column]] %in% columns_unselected,])
     if (length(input$selected_samples) > 0){
      data$retain_samples <- retain_samples
      updateSelectInput(session, "remove_samples", choices=c(retain_samples), selected=input$remove_samples)
      updateSelectInput(session, "include_samples", choices=c(rest_of_samples), selected=input$include_samples)
    } else {
      data$retain_samples  <- NULL
      updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=input$include_samples)

    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)


  observeEvent(input$remove_samples, {
    if ( !is.null(input$remove_samples)){
      outliersamples <- input$remove_samples
      data$remove_samples <- outliersamples
    } else {
      data$remove_samples  <- NULL
    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)


  observeEvent(input$include_samples, {
    if ( length(input$include_samples) > 0){
      includesamples <- input$include_samples
      data$include_samples <- includesamples
    } else {
      data$include_samples  <- NULL
    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)


   timer <- reactiveValues(data = NULL, redraw = TRUE)


    observe({
      data$retain_samples
      data$include_samples
      data$remove_samples
      timer$redraw <- FALSE
    })


    observe({
      invalidateLater(1000, session)
      data$retain_samples
      data$include_samples
      data$remove_samples
      if (isolate(timer$redraw)) {
        timer$data <- subset(fusions$summary, sample %in% setdiff(unique(c(data$retain_samples, data$include_samples)),data$remove_samples))
      } else {
        isolate(timer$redraw <- TRUE)
      }
    })


    observe({
      if( !(length(input$selected_samples) > 0) & !( length(input$include_samples) > 0) ){
        timer$data <- NULL
      }
    })



   ## plot --------------------------------------------
   output$overview <- renderPlot({

     shiny::validate(
       need(nrow(timer$data) > 0,
            message = "Select samples to get started.")
     )

     if ( length(unique(timer$data$variable)) > length(ggsci_db[[input$color_palette]]$default) ){
       color_set = colorRampPalette(ggsci_db[[input$color_palette]]$default)( length(unique(timer$data$variable)) )
     } else {
       color_set = as.vector( ggsci_db[[input$color_palette]]$default[1:length(unique(timer$data$variable)) ])
     }

     P <- ggplot( timer$data ,aes(x = sample, y = value,fill = variable)) + theme_bw(base_size=10) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  coord_flip() + xlab("") + ylab("") + ggtitle(input$title) + theme(plot.title = element_text(size=16), legend.position=input$legend_position)

     if (input$percentage){
       P <- P +   geom_bar(position = "fill",stat = "identity") +   scale_y_continuous(labels = percent_format())
     } else {
      P <- P +   geom_bar(stat = "identity")
     }


     if(input$color_palette != "ggplot"){
       P <- P + scale_fill_manual(values=color_set)
     }

     facets <- paste(input$facet_by_row, '~', input$facet_by_col)
     if (facets != '. ~ .')
       P <- P + facet_grid(facets, scales="free")



      if(!input$sample_labels)
      P <- P +  theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank())

     P

   })


   observe({
     updateSelectizeInput(session, "selected_fusion", choices = c(unique(subset(fusions$list, sample %in% setdiff(unique(c(data$retain_samples, data$include_samples)),data$remove_samples))$fusion.name)),selected=NULL, server=TRUE)
   })




    output$fusion_table <-  DT::renderDataTable({

      df <- subset(fusions$list, sample %in% setdiff(unique(c(data$retain_samples, data$include_samples)),data$remove_samples))
      if(length(input$selected_fusion) > 0 ){
        df <- subset(df, fusion.name  %in% c(input$selected_fusion))
      }
      as.data.frame(df)# [,-c(3,10,11,12,13)]
    },server=TRUE, rownames = FALSE, filter = list(position = 'top', clear = FALSE), options = list(scrollX = TRUE,
                     search=list(regex = TRUE, caseInsensitive = TRUE),
                     pageLength = 20             ) )
# 'fusion.name' = 9,




    observeEvent({

      shiny::validate(
        need(input$tabs == "Intersections",
             message = "")
      )


      input$select_samples_by
      input$set_size    }, {
        if (input$select_samples_by != ""){

  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, {
                 for (i in 1:4) {
                   list_of_fusions = lapply(unique(fusions$list[[input$select_samples_by]]),function(x)  fusions$list[get(input$select_samples_by) == x]$fusion.name )
                   incProgress(2/4)
                   names(list_of_fusions) <- unique(fusions$list[[input$select_samples_by]])
                   incProgress(1/4)
                   vennset <- systemPipeR::overLapper(list_of_fusions, type="vennsets")
                   incProgress(1/4)
                   Sys.sleep(0.25)
                 }
               })

  l <- systemPipeR::vennlist(vennset)
  l <- l[ lengths(l) > input$set_size ]
  max_table = length(l)
    # https://stackoverflow.com/questions/28003715/caption-in-rendertable-shiny

    lst <- list()
    for (i in 1:max_table) {
      df <- subset(fusions$list, fusion.name %in% l[[i]])
      lst[[i]] <-  as.data.frame(df)#[,-c(3,10,11,12,13)]

    }


    output$inter_tables <- renderUI({
      plot_output_list <- lapply(1:max_table, function(i) {
        tablename <- paste("tablename", i, sep = "")
        tabPanel(names(l[i]), DT::dataTableOutput(session$ns(tablename)))
      })
      do.call(tabsetPanel,c(plot_output_list,id="level"))
    })


    for (i in 1:max_table) {
      local({
        my_i <- i
        tablename <- paste("tablename", my_i, sep = "")
        output[[tablename]] <- DT::renderDataTable({
          lst[[my_i]]
        },server=TRUE, rownames = FALSE, filter = list(position = 'top', clear = FALSE), options = list(scrollX = TRUE,
                                                                                                                                                                                                                                                                             search=list(regex = TRUE, caseInsensitive = TRUE),
                                                                                                                                                                                                                                                                             pageLength = 20             ) )
      })
    }
}
    }, ignoreNULL =  FALSE)




    output$upsetr <- renderPlot({


    })







    # plot --------------------------------------------
    output$circos_plot <- renderPlot({

      shiny::validate(
        need(nrow(timer$data) > 0,
             message = "Select samples to get started.")
      )


      withProgress(message = 'Processing',
                   detail = 'This may take a moment...', value = 0, {

                     fus_sub <- as.data.frame(subset(fusions$list, sample %in% setdiff(unique(c(data$retain_samples, data$include_samples)),data$remove_samples)))
                     #  fus_sub <- as.data.frame(head(fusions$list, 20))
                     which_method = input$use_fusions_from # "star.fusion_final"
                     fus_sub <- fus_sub[, c("Left.genes", "Right.genes", which_method)]

                     fus_sub <- fus_sub %>%  dplyr::mutate(fus = strsplit(as.character(as.character(fus_sub[,grep(which_method, colnames(fus_sub))])), ",")) %>% tidyr::unnest(fus)
                     fus_sub <- fus_sub[,-3]
                     fus_sub$reads <-  stringr::str_split_fixed(fus_sub$fus, "=", 2)[,1]
                     fus_sub$chromLeftStart <- as.matrix(stringr::str_split_fixed(fus_sub$fus, "=", 2)[,2] %>% stringr::str_split_fixed(., "-", 2))[,1]
                     fus_sub$chromLeft <- stringr::str_split_fixed(fus_sub$chromLeftStart, ":", 2)[,1]
                     fus_sub$chromRightStart <- as.matrix(stringr::str_split_fixed(fus_sub$fus, "=", 2)[,2] %>% stringr::str_split_fixed(., "-", 2))[,2]
                     fus_sub$chromRight <- stringr::str_split_fixed(fus_sub$chromRightStart, ":", 2)[,1]
                     fus_sub$chromLeftStart <- gsub(".*\\:","",fus_sub$chromLeftStart)
                     fus_sub$chromRightStart <- gsub(".*\\:","",fus_sub$chromRightStart)
                     fus_sub <- fus_sub[,-3]

                     geneLabelDataLeft <- data.frame(Chromosome=fus_sub$chromLeft, chromStart=as.numeric(fus_sub$chromLeftStart), chromEnd=as.numeric(fus_sub$chromLeftStart)+1, Gene=fus_sub$Left.genes)
                     geneLabelDataRight <- data.frame(Chromosome=fus_sub$chromRight, chromStart=as.numeric(fus_sub$chromRightStart), chromEnd=as.numeric(fus_sub$chromRightStart)+1, Gene=fus_sub$Right.genes)
                     geneLabelData <-  unique(rbind(geneLabelDataLeft, geneLabelDataRight))

                     cytobandFile <- system.file("extdata", "UCSC.HG19.Human.CytoBandIdeogram.txt",  package = "chimeraviz")
                     cytoband <- utils::read.table(cytobandFile)
                     names(cytoband) <- c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain")
                     assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)
                     cytoband <- RCircos::RCircos.Sort.Genomic.Data(genomic.data = cytoband, is.ideo = TRUE)
                     cyto.info <- cytoband
                     chr.exclude <- NULL
                     tracks.inside <- 3
                     tracks.outside <- 0
                     RCircos::RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
                     params <- RCircos::RCircos.Get.Plot.Parameters();
                     params$text.size <- 0.8
                     params$char.width <- 500
                     params$text.color <- "black"
                #     params$margins = 1
                     RCircos::RCircos.Reset.Plot.Parameters(params);
                     RCircos.Env$RCircos.PlotPar$text.size = .8
                     RCircos::RCircos.Set.Plot.Area(margins = 0)
                     RCircos::RCircos.Chromosome.Ideogram.Plot()

                     name.col <- 4
                     side <- "in"
                     track.num <- 1
                     RCircos::RCircos.Gene.Connector.Plot(geneLabelData, track.num, side)
                     track.num <- 2
                     RCircos::RCircos.Gene.Name.Plot(geneLabelData, name.col, track.num, side)

                     linkData <-  data.frame(Chromosome=fus_sub$chromLeft, chromStart=as.numeric(fus_sub$chromLeftStart), chromEnd=as.numeric(fus_sub$chromLeftStart)+1, Chromosome.1=fus_sub$chromRight, chromStart.1=as.numeric(fus_sub$chromRightStart), chromEnd.1=as.numeric(fus_sub$chromRightStart)+1,  linkWidth=as.numeric(fus_sub$reads))
                     linkData$linkWidth <- chimeraviz:::.scale_list_to_interval(linkData$linkWidth, 1,6)

                     track.num <- 3
                     RCircos::RCircos.Link.Plot(link.data = linkData, track.num = track.num,
                                                by.chromosome = TRUE, start.pos = NULL, genomic.columns = 3,
                                                is.sorted = FALSE, lineWidth = linkData$linkWidth)

                     title(input$title, line = -5)
                     remove("RCircos.Env", envir = .GlobalEnv)

                   })

    })









}