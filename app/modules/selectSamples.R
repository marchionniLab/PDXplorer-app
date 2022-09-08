

selectSamples <- function(input, output, session, dds, metadata) {

  data <- reactiveValues(dds = NULL, retain_samples = NULL, remove_samples = NULL, include_samples = NULL, redraw=FALSE, genedata= NULL, final_list_samples_to_include = NULL)

  updateSelectInput(session, "remove_samples", choices=c("",colnames(dds)), selected="")
  updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected="type")
  updateSelectInput(session, "columun_annotation", choices = c("", names(metadata)), selected="")
  updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=NULL)

  observeEvent(input$select_samples_by, {
    updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected=input$selected_samples, inline=TRUE )
  }, ignoreInit = TRUE)


  observeEvent(input$select_all,{
    updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= unique(metadata[[input$select_samples_by]]), inline=TRUE )
  })



  observeEvent(input$selected_samples, {
    if ( length(input$selected_samples) > 0){
      selected_column <- input$select_samples_by
      selected_groups <- input$selected_samples
      retain_samples <- row.names(metadata[metadata[[selected_column]] %in% selected_groups,])
      columns_unselected <- setdiff(unique(metadata[[selected_column]]) , selected_groups)
      rest_of_samples <- row.names( metadata[metadata[[selected_column]] %in% columns_unselected,])
      updateSelectInput(session, "remove_samples", choices=c(retain_samples), selected=input$remove_samples)
      updateSelectInput(session, "include_samples", choices=c(rest_of_samples), selected=input$include_samples)
      data$retain_samples <- retain_samples
    } else {
      updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=input$include_samples)
      data$retain_samples <- NULL
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
    if ( !is.null(input$include_samples)){
      includesamples <- input$include_samples
      data$include_samples <- includesamples
    } else {
      data$include_samples  <- NULL
    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)


  observe({
    data$final_list_samples_to_include <- setdiff(colnames(dds[,unique(c(data$retain_samples, data$include_samples))]),data$remove_samples)
  })


  timer <- reactiveValues(final_list_samples_to_include=NULL)


  observe({
    data$final_list_samples_to_include
    timer$redraw <- FALSE
  })


  observe({
    invalidateLater(1000, session)
    data$final_list_samples_to_include
    if (isolate(timer$redraw)) {
       timer$final_list_samples_to_include <- data$final_list_samples_to_include
    } else {
      isolate(timer$redraw <- TRUE)
    }
  })


}
