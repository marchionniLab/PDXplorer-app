
geInput <- function(id) {
  ## Define what's shown on the side for all plots of the Gene Exploration panel

  ns <- NS(id) # namespaced IDs; here: "ge" is a likely candidate

  title <- textInput(ns("title"),
                     label = "Plot title",
                     value = "")

  gene_symbol <- selectizeInput(ns("gene_symbol"),
                              "Select gene to plot",
                              choices = NULL,
                              selected=NULL,
                              multiple=TRUE)

  paste_from_clipboard <- textAreaInput(ns('paste_from_clipboard'), label = 'Include below genes: ', value = "")


  multi_gene_symbol <- selectizeInput(ns("multi_gene_symbol"),
                                "Select genes to plot",
                                choices = NULL,
                                selected=NULL,
                                multiple=TRUE)
  top_n_genes <- numericInput(ns('top_n_genes'), label = 'Include top n most variable genes: ', value = 0, min = 0,max = 5000, step=25)

  dist_method <- selectInput(inputId=ns("dist_method"),
                             label="Select distance measure",
                             choices = c("euclidean",  "maximum", "manhattan", "canberra", "minkowski"),
                             selected = "euclidean")

  plot_style <- selectInput(ns('plot_style'),
                            "Plot style for gene count",
                            choices = list("boxplot","violin plot"))

  select_all <- actionLink(ns("select_all"), "Select all", style = "color:#18cc7b")
  deselect_all <- actionLink(ns("deselect_all"), "Deselect all", style = "color:#cc3c18")

  group_by <- selectizeInput(inputId=ns("group_by"),
                             label="Group by",
                             choices = NULL,
                             selected=NULL,
                             multiple=TRUE)


  color_by <- selectizeInput(inputId=ns("color_by"),
                             label="Color by",
                             choices = NULL,
                             selected=NULL,
                             multiple=FALSE)

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

  remove_batch_effect <-selectInput(ns('remove_batch_effect'), label = "Remove batch effect?", choices = NULL, selected=NULL, multiple = FALSE)


  remove_samples <-selectInput(ns('remove_samples'), label = "Remove sample by name:", choices = NULL, selected=NULL, multiple = TRUE)

  include_samples <-selectInput(ns('include_samples'), label = "Include sample by name:", choices = NULL, selected=NULL, multiple = TRUE)


  # list of UIs to be exported to the add
  tagList(title,
          gene_symbol,
          multi_gene_symbol,
          paste_from_clipboard,
       top_n_genes,
          dist_method,
          group_by,
       color_by,
          plot_style,
          facet_by_row,
          facet_by_col,
       hr(),
       remove_batch_effect,
          hr(),
          h5("Include / exclude samples"),
          select_samples_by,
       tags$div(select_all, style = "display:inline-block; float:left"),
       tags$div(deselect_all,  style = "display:inline-block; float:right"),
       br(),
          selected_samples,
          remove_samples,
          include_samples)
 #         facet_by_col

}


geOutput <- function(id) {
  ## Define the different types of plots that will be shown
  ns <- NS(id)


  tagList(
    ## cluster
    tabsetPanel(id = ns("tabs"),
                tabPanel("Boxplots",

    shinydashboard::box( width = 12, title = "Gene expression boxplot",
                         fluidRow(column(width = 3,
                                         checkboxInput(ns('sample_labels'),
                                                       "Display sample labels?",
                                                       value = FALSE)),
                                  column(width=3,
                                         checkboxInput(ns('ylimZero'),
                                                       "Set y axis limit to 0?",
                                                       value = TRUE)),
                                  column(width=3,
                                         checkboxInput(ns('rotate_axis'),
                                                       "Rotate X-axis labels?",
                                                       value = TRUE)),
                                  column(width=3,
                                         checkboxInput(ns('show_dots'),
                                                       "Show dots?",
                                                       value = FALSE))),
                         fluidRow(
                           column(width = 3,
                                  selectInput(ns('color_palette'),
                                              label = 'Select color palette:',
                                              choices = c(names(ggsci_db)[-1]),
                                              selected = "npg" ,multiple = FALSE)),
                           column(width = 3,
                                selectInput(ns('legend_position'), label = 'Legend position:', choices = c("bottom", "top", "right", "left", "none" = "none"),
                                            selected = "right",multiple = FALSE))),

                         plotOutput(ns("gene_symbol_plot"), height="800px", click = clickOpts(id =ns("plot_click_ge"))),
                         verbatimTextOutput(ns("click_info")),
						 ggDownloadUI(ns("ggbox")),
                         collapsible = FALSE)
                ),
    tabPanel("Heatmap",
             shinydashboard::box( width = 12, title = "Gene expression heatmap",
                                  fluidRow(column(width = 3,
                                                  checkboxInput(ns('sample_labels_hm'),
                                                                "Display sample labels?",
                                                                value = TRUE)),
                                           column(width = 3,
                                                  checkboxInput(ns('row_labels_hm'),
                                                                "Display gene IDs?",
                                                                value = TRUE)),
                                           column(width=3,
                                                  checkboxInput(ns('cluster_col'),
                                                                "Cluster columns?",
                                                                value = TRUE)),
                                           column(width=3,
                                                  checkboxInput(ns('cluster_row'),
                                                                "Cluster rows?",
                                                                value = TRUE))

                                           ),
                                  fluidRow(column(width = 3,
                                                  selectInput(ns('column_annotation'), label = 'Annotate column by:', choices = NULL,
                                                              selected = NULL,multiple = TRUE)),
                                           column(width = 3,
                                                  selectInput(ns('cutree_cols'), label = 'Cut cols into N groups:', choices = seq(1,20), selected = 1,
                                                              multiple = FALSE)),
                                           column(width = 3,
                                                  selectInput(ns('scale_by'),
                                                              label = 'Scale by:',
                                                              choices = c("none", "row", "column"),
                                                              selected = "row" ,multiple = FALSE)),
                                           column(width = 3,
                                                  selectInput(ns('color_palette_seq'),
                                                              label = 'Select color palette:',
                                                              choices = row.names(subset(RColorBrewer:::brewer.pal.info, category=="div"  | category=="seq")),
                                                              selected = "RdYlBu" ,multiple = FALSE))
                                  ),
                                  plotOutput(ns("pheatmap_plot"), height="800px")),
								  ggDownloadUI(ns("ggheat")),
             collapsible = FALSE
             )

    )

  )

}


geMod <- function(input, output, session, dds = dds, metadata) {

  updateSelectizeInput(session, "gene_symbol", choices = c(row.names(dds)), selected=row.names(dds[1,]), server=TRUE)
  updateSelectizeInput(session, "multi_gene_symbol", choices = c(row.names(dds)), selected=NULL, server=TRUE)
  updateSelectInput(session, "group_by", choices = c("", names(metadata), "gene"), selected = "type")
  updateSelectInput(session, "color_by", choices = c("", names(metadata), "gene"), selected = "type")
  updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected="type")
  updateSelectInput(session, "facet_by_row", choices = c(none = ".", names(metadata), "gene"), selected=".")
  updateSelectInput(session, "facet_by_col", choices = c(none = ".", names(metadata), "gene"), selected=".")
  updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=NULL)
  updateSelectInput(session, "column_annotation", choices = c("", names(metadata)), selected="type")
  updateSelectInput(session, "remove_batch_effect", choices = c("none","source", "batch", "patient", "purity"), selected=NULL)



  observeEvent(input$tabs, {
    if (input$tabs ==  "Boxplots"){
      shinyjs::hide("multi_gene_symbol")
      shinyjs::show("gene_symbol")
      shinyjs::show("gene_symbol")
      shinyjs::show("group_by")
      shinyjs::show("color_by")
      shinyjs::show("plot_style")
      shinyjs::show("facet_by_row")
      shinyjs::show("facet_by_col")
      shinyjs::hide("dist_method")
      shinyjs::hide("remove_batch_effect")
      shinyjs::hide("paste_from_clipboard")
      shinyjs::hide("top_n_genes")
    } else if (input$tabs ==  "Heatmap"){
      shinyjs::show("top_n_genes")
      shinyjs::show("remove_batch_effect")
      shinyjs::show("multi_gene_symbol")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("group_by")
      shinyjs::hide("plot_style")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("show")
      shinyjs::hide("color_by")
      shinyjs::show("paste_from_clipboard")
    }
  })

  mydds <- reactiveValues(retain_samples = NULL, remove_samples = NULL, include_samples = NULL, final_list_samples_to_include=NULL, from_paste=NULL)
  timer <- reactive(list(final_list_samples_to_include = mydds$final_list_samples_to_include,
                         retain_samples  = mydds$retain_samples,
                         include_samples = mydds$include_samples,
                         from_paste=mydds$from_paste,
                         multi_gene_symbol=input$multi_gene_symbol,
                         top_n_genes=input$top_n_genes)) %>% debounce(1000)

  observeEvent(input$select_samples_by, {
    updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= unique(metadata[[input$select_samples_by]]), inline=TRUE )

  }, ignoreInit = FALSE)
  observeEvent(input$select_all,{
    updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= unique(metadata[[input$select_samples_by]]), inline=TRUE )
  })
  observeEvent(input$deselect_all,{
    updateCheckboxGroupInput(session, "selected_samples",  choices =  unique(metadata[[input$select_samples_by]]), selected= NULL, inline=TRUE )
  })
  observeEvent(input$selected_samples, {
    selected_column <- input$select_samples_by
    selected_groups <- input$selected_samples
    retain_samples <- row.names(metadata[metadata[[selected_column]] %in% selected_groups,])

    columns_unselected <- setdiff(unique(metadata[[selected_column]]) , selected_groups)
    rest_of_samples <- row.names( metadata[metadata[[selected_column]] %in% columns_unselected,])
    if (length(input$selected_samples) > 0){
      updateSelectInput(session, "remove_samples", choices=c(retain_samples), selected=input$remove_samples)
      updateSelectInput(session, "include_samples", choices=c(rest_of_samples), selected=input$include_samples)
      mydds$retain_samples <- retain_samples
    } else {
      mydds$retain_samples <- NULL
      updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=input$include_samples)

    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)
  observeEvent(input$remove_samples, {
    if ( !is.null(input$remove_samples)){
      outliersamples <- input$remove_samples
      mydds$remove_samples <- outliersamples
    } else {
      mydds$remove_samples  <- NULL
    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)
  observeEvent(input$include_samples, {
    if ( !is.null(input$include_samples)){
      includesamples <- input$include_samples
      mydds$include_samples <- includesamples
    } else {
      mydds$include_samples  <- NULL
    }
  }, ignoreInit = FALSE, ignoreNULL=FALSE)
  observe({
    mydds$final_list_samples_to_include <- setdiff(colnames(dds[,unique(c(mydds$retain_samples, mydds$include_samples))]),mydds$remove_samples)
  })

  gene_data <- reactive({
    shiny::validate(
      need(length(timer()$retain_samples) >0  | length(timer()$include_samples) > 0,
           message = "Select a sample to get started.")
    )

    shiny::validate(
      need(input$gene_symbol != "",
           message = "Select a gene to get started.")
    )

    shiny::validate(
      need(input$group_by != "",
           message = "Select a factor to group samples by.")
    )


    withProgress(message = 'Extracting gene data', {
      incProgress(0)
      object = dds[,timer()$final_list_samples_to_include]
   #  genedata <- plotCounts(object, gene=input$gene_symbol,intgroup = colnames(colData(object)),returnData = TRUE)


      genedata <- assay(object[input$gene_symbol,]) %>% reshape2::melt(., id.vars = rownames)
      names(genedata) <- c("gene","sample","count")
      genedata <- merge(genedata, colData(object), by.x="sample", by.y="sampleName")
      genedata <- as.data.frame(genedata)


      return(genedata)
    })
  }) %>%
    bindCache(timer()$final_list_samples_to_include,input$gene_symbol, input$group_by)



	makeGenePlot  <- function() {
	    shiny::validate(
	      need( !is.null(gene_data()),
	            message = "Select samples to get started.")
	    )

	    shiny::validate(
	      need(length(timer()$final_list_samples_to_include) > 0,
	           message = "Select samples to get started.")
	    )

	    shiny::validate(
	      need(input$gene_symbol != "",
	           message = "Select a gene to get started.")
	    )

	    shiny::validate(
	      need(input$group_by != "",
	           message = "Select a factor to group samples by.")
	    )




	    genedata = gene_data()
	    factors <- genedata[,match(input$group_by,colnames(genedata))]
	    genedata$group <- interaction(factors)

	    if(input$color_by != "none"){
	      if ( length(unique(genedata[[input$color_by]])) > length(ggsci_db[[input$color_palette]]$default) ){
	        color_set = colorRampPalette(ggsci_db[[input$color_palette]]$default)(  length(unique(genedata[[input$color_by]]))  )
	      } else {
	        color_set = as.vector( ggsci_db[[input$color_palette]]$default[1: length(unique(genedata[[input$color_by]]))])
	      }
	    }



	    if(input$plot_style=="boxplot"){
	      P <- ggplot(genedata,aes_string(x="group",y="count",fill=input$color_by, linetype="gene")) +
	        geom_boxplot(outlier.shape = NA,alpha=0.7) + theme_bw(base_size=14)
	      if(input$ylimZero){
	        P <- P + scale_y_log10(name="Normalized counts - log10 scale", limits=c(0.4,NA))
	      } else {
	        P <- P + scale_y_log10(name="Normalized counts - log10 scale")
	      }

	      if(input$show_dots==TRUE){
	        P <- P +
	          labs(title=paste0("Normalized counts for "," - ",paste(input$gene_symbol, collapse=","))) +
	          scale_x_discrete(name="") +
	          geom_jitter(aes_string(x="group",y="count"),position = position_jitter(width = 0.1)) +
	          scale_fill_discrete(name="Experimental\nconditions")
	      } else {
	        P <- P +
	          labs(title=paste0("Normalized counts for "," - ",paste(input$gene_symbol, collapse=","))) +
	          scale_x_discrete(name="") +
	          scale_fill_discrete(name="Experimental\nconditions")
	      }
	      if(input$sample_labels == TRUE){
	        P <- P + geom_text(aes_string(label="sampleName"),hjust=-.1,vjust=0)
	      }
	    } else if(input$plot_style=="violin plot"){
	      P <- ggplot(genedata,aes_string(x="group",y="count",fill=input$color_by, linetype="gene")) +
	        geom_violin(aes_string(col="group"),alpha = 0.6) +  theme_bw(base_size=14)
	      if(input$ylimZero){
	        P <- P + scale_y_log10(name="Normalized counts - log10 scale", limits=c(0.4,NA))
	      } else {
	        P <- P + scale_y_log10(name="Normalized counts - log10 scale")
	      }

	      if(input$show_dots==TRUE){
	        P <- P +
	          labs(title=paste0("Normalized counts for "," - ",paste(input$gene_symbol, collapse=","))) +
	          scale_x_discrete(name="") +
	          geom_jitter(aes_string(x="group",y="count"),alpha = 0.8,position = position_jitter(width = 0.1)) +
	          scale_fill_discrete(name="Experimental\nconditions") + scale_color_discrete(guide="none")  + scale_fill_manual(values=color_set)
	      } else {
	        P <- P +
	          labs(title=paste0("Normalized counts for "," - ",paste(input$gene_symbol, collapse=","))) +
	          scale_x_discrete(name="") +
	          scale_fill_discrete(name="Experimental\nconditions") + scale_color_discrete(guide="none")  + scale_fill_manual(values=color_set)
	      }

	      if(input$sample_labels){
	        P <- P + geom_text(aes_string(label="sampleName"),hjust=-.1,vjust=0)
	      }
	    }

	    if(input$rotate_axis){
	      P <- P + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	    }
	    facets <- paste(input$facet_by_row, '~', input$facet_by_col)
	    if (facets != '. ~ .')
	      P <- P + facet_grid(facets, scales="free")

	     P <- P + ggtitle(input$title)  +  theme(legend.position=input$legend_position)

	    if(input$color_palette != "ggplot"){
	      P <- P + scale_fill_manual(values=color_set)
	    }

	     if(input$title == "")
	       P <- P + ggtitle(paste0("Normalized counts for "," - ",paste(input$gene_symbol, collapse=",")))
	     P

	}

  output$gene_symbol_plot <- renderCachedPlot({

    shiny::validate(
      need( !is.null(gene_data()),
            message = "Select samples to get started.")
    )

    shiny::validate(
      need(length(timer()$final_list_samples_to_include) > 0,
           message = "Select samples to get started.")
    )

    shiny::validate(
      need(input$gene_symbol != "",
           message = "Select a gene to get started.")
    )

    shiny::validate(
      need(input$group_by != "",
           message = "Select a factor to group samples by.")
    )


	  makeGenePlot()

  }, cacheKeyExpr = { list(gene_data(),
              input$gene_symbol,
              input$title,
              input$legend_position,
              input$facet_by_row,
              input$facet_by_col,
              input$rotate_axis,
              input$show_dots,
              input$sample_labels,
              input$plot_style,
              input$ylimZero,
              timer()$final_list_samples_to_include,
              input$color_palette,
              input$color_by,
              input$group_by)})





  observeEvent(input$paste_from_clipboard,{
    from_clip <- input$paste_from_clipboard
    input_genes = c(trimws(as.character(strsplit(from_clip, "\n")[[1]])))
    input_genes = stringr::str_to_upper(unique(input_genes))
    input_genes = input_genes[input_genes %in% row.names(dds)]
    mydds$from_paste <-  input_genes
  })



  heatmap_data <- reactive({
    shiny::validate(
      need(length(timer()$final_list_samples_to_include) >0,
           message = "Select a sample to get started.")
    )

    shiny::validate(
      need(timer()$multi_gene_symbol != "" || timer()$top_n_genes > 0 || timer()$from_paste != "",
           message = "Select a gene to get started.")
    )


    withProgress(message = 'Extracting heatmap data', {
      incProgress(0)
      object = dds[,timer()$final_list_samples_to_include]

      vst_data <- assay(object, "vst")

      if (input$remove_batch_effect != "none"){
        tryCatch(vst_data <- limma::removeBatchEffect(vst_data, batch=colData(object)[[input$remove_batch_effect]]), error=function(e) vst_data <- assay(object, "vst"))
      }

      if(length(timer()$from_paste) > 0){
        from_paste = timer()$from_paste
      } else {
        from_paste <- NULL
      }


      if(timer()$top_n_genes > 0){
        ntop = timer()$top_n_genes
        rv <- rowVars(vst_data)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        genes_from_var <- row.names(vst_data[select,])
      } else {
        genes_from_var <- NULL
      }

      gene_symbols <- unique(c(timer()$multi_gene_symbol, genes_from_var, from_paste))

      shiny::validate(
        need(length(gene_symbols) > 1,
             message = "Select at least two genes to get started.")
      )


      pheatmap_data <- vst_data[gene_symbols, ]

      return(pheatmap_data)
    })
  }) %>%
    bindCache(timer()$final_list_samples_to_include,timer()$multi_gene_symbol, timer()$top_n_genes, input$remove_batch_effect, timer()$from_paste)




		makeHeatPlot  <- function() {

		    shiny::validate(
		      need(length(timer()$final_list_samples_to_include) > 0,
		           message = "Select samples to get started.")
		    )

		    shiny::validate(
		      need(nrow(heatmap_data()) > 1,
		           message = "Select at least two genes to get started.")
		    )


		    heatmap_data <- heatmap_data()


		    if(input$scale_by !="none"){
		      heatmap_data <- heatmap_data[apply(heatmap_data, MARGIN = 1, FUN = function(x) sd(x) != 0),]
		      shiny::validate(
		        need( nrow(heatmap_data) > 1,
		              message = paste("Turn off row scaling to view selected genes.", "The following genes were selected:", paste(row.names(heatmap_data), collapse=", ") ))
		      )
		    }


		    if (length(input$column_annotation) > 0){
		      expgroups <- as.data.frame(colData(dds)[,input$column_annotation])
		      rownames(expgroups) <- colnames(dds)
		      colnames(expgroups) <- input$column_annotation
		    } else {
		      expgroups = NULL
		    }


		    pheatmap::pheatmap(heatmap_data,
		                       cluster_rows = input$cluster_row,
		                       cluster_cols = input$cluster_col,
		                       clustering_distance_rows=input$dist_method,
		                       clustering_distance_cols=input$dist_method,
		                       show_rownames = input$row_labels_hm,
		                       show_colnames = input$sample_labels_hm,
		                       cutree_cols = input$cutree_cols,
		                       main=input$title,
		                       annotation_col = expgroups,
		                       scale=input$scale_by,
		                       color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100))

		}

  output$pheatmap_plot <- renderCachedPlot({

    shiny::validate(
      need(length(timer()$final_list_samples_to_include) > 0,
           message = "Select samples to get started.")
    )

    shiny::validate(
      need(nrow(heatmap_data()) > 1,
           message = "Select at least two genes to get started.")
    )

makeHeatPlot()

  }, cacheKeyExpr = { list(heatmap_data(),
                           timer()$final_list_samples_to_include,
                           timer()$top_n_genes,
                           timer()$multi_gene_symbol,
                           timer()$from_paste,
                           input$remove_batch_effect,
                           input$color_palette_seq,
                           input$scale_by,
                           input$cutree_cols,
                           input$sample_labels_hm,
                           input$row_labels_hm,
                           input$title,
                           input$dist_method,
                           input$cluster_row,
                           input$cluster_col,
                           input$column_annotation,
                           input$scale_by)})


callModule(ggDownload, "ggbox", reactive({ makeGenePlot()}))


callModule(ggDownload, "ggheat", reactive({makeHeatPlot()}))

}




