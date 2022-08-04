markersInput <- function(id) {

  ns <- NS(id)


#  subset_by <- selectInput(inputId=ns("subset_by"),
 #                                  label="Subset by",
#                                   choices = c("none", "PDX", "primary"),
#                                   selected = "none")
  title <- textInput(ns("title"), label = "Plot title", value = "")


  find_markers_by <- selectInput(inputId=ns("find_markers_by"),
                                   label="Find markers by",
#                                   choices = c("type", "source", "PDX : type", "primary : type"),
                             choices = c("type"),
                                   selected = "type")


  direction <- selectInput(ns('direction'), label = "Direction of log-fold changes to be considered:", choices=c("any", "up", "down"), selected="any")


  top_genes <- numericInput(ns('top_genes'), label = "Use genes in the top n rank:", min = 1, max = 25, value = 5)

  pvalue_type <- numericInput(ns('pvalue_type'), label = "Use genes in the top n rank:", min = 1, max = 10, value = 5)

  remove_duplicates <- checkboxInput(ns('remove_duplicates'), "Remove duplicate genes?", value = FALSE)

  tagList(#subset_by,
    title,
          find_markers_by,
          hr(),
          h5("Filter results"),
          direction,
          top_genes#,
        #  remove_duplicates
  )

}



markersOutput <- function(id) {
  ns <- NS(id)

  tagList(

    ## cluster
    tabsetPanel(id = ns("meta_tabs"),
                tabPanel("Marker genes",
                shinydashboard::box( width = 12, title = "Marker genes",
                                     fluidRow(column(width = 2,
                                                     checkboxInput(ns('sample_labels_hm'),
                                                                   "Display sample labels?",
                                                                   value = TRUE)),
                                              column(width = 2,
                                                     checkboxInput(ns('row_labels_hm'),
                                                                   "Display gene IDs?",
                                                                   value = TRUE)),
                                              column(width=2,
                                                     checkboxInput(ns('cluster_col'),
                                                                   "Cluster columns?",
                                                                   value = TRUE)),
                                              column(width=2,
                                                     checkboxInput(ns('cluster_row'),
                                                                   "Cluster rows?",
                                                                   value = TRUE)),
																   column(width=2,
								                                          checkboxInput(ns('remove_source_effect'),
								                                                        "Remove source effect?",
								                                                        value = FALSE))

                                     ),
                                     fluidRow(
                                   column(width = 3,
                                                       selectInput(ns('column_annotation'), label = 'Annotate column by:', choices = NULL,
                                                                   selected = NULL,multiple = TRUE)),
						                                           column(width = 2,
						                                                  selectInput(ns('scale_by'),
						                                                              label = 'Scale by:',
						                                                              choices = c("none", "row", "column"),
						                                                              selected = "row" ,multiple = FALSE)),
						                                           column(width = 2,
						                                                  selectInput(ns('color_palette_seq'),
						                                                              label = 'Select color palette:',
						                                                              choices = row.names(subset(RColorBrewer:::brewer.pal.info, category=="div"  | category=="seq")),
						                                                              selected = "RdYlBu" ,multiple = FALSE)),
                                        column(width = 2,
                                                      numericInput(ns('fontsize_row'), label = 'Row fontsize :', value = 8,
                                                                  min = 1,step = 1)
                                                      ),

                                              column(width = 2,
                                                     selectInput(ns('cutree_cols'), label = 'Cut cols into N groups:', choices = seq(1,20), selected = 1,
                                                                 multiple = FALSE))
                                     ),
                                     plotOutput(ns("pheatmap_plot"), height="800px"),
                                     textOutput(ns('gene_list')),
                                     helpText("Using a modified findMarkers function from the scran package,
this function uses limma to test for differentially expressed genes (DEGs) between groups of samples. For each group, the log-fold changes and other statistics from all relevant pairwise comparisons are combined into a single table. A list of such tables is returned for all groups to define a set of potential marker genes.  Each table is sorted by the Top value, which specifies the size of the candidate marker set. Taking all rows with Top values no greater than some integer X will yield a set containing the top X genes (ranked by significance) from each pairwise comparison. For example, if X is 5, the set will consist of the union of the top 5 genes from each pairwise comparison. The marker set for each cluster allows it to be distinguished from the other clusters based on the expression of at least one gene.")),
ggDownloadUI(ns("ggheat")),
                collapsible = FALSE
    )
    )
  )
}


markersMod <- function(input, output, session, dds, metadata) {

  updateSelectInput(session, "column_annotation", choices = c("", names(metadata)), selected="type")


  markers_data <- reactive({
      shiny::validate(
        need(input$find_markers_by != "",
             message = "")
      )

    withProgress(message = 'Extracting heatmap data', {
      incProgress(0)
      if (input$find_markers_by == "source"){
        counts <- assay(dds, "counts")
        selected_meta <- metadata
        retain <- row.names(subset(metadata, source == "PDX" | source == "primary"))
        clusters <- factor(selected_meta$source)
        design=NULL
        column_annotation = "source"
      } else  if (input$find_markers_by == "PDX : type"){
        retain <- row.names(subset(metadata, source == "PDX"))
        counts <- assay(dds[,retain], "counts")
        selected_meta <- metadata[retain,]
        clusters <- factor(selected_meta$type)
        design=NULL
        column_annotation = "type"
      } else if (input$find_markers_by == "primary : type"){
        retain <- row.names(subset(metadata, source == "primary"))
        counts <- assay(dds[,retain], "counts")
        selected_meta <- metadata[retain,]
        clusters <- factor(selected_meta$type)
        column_annotation = "type"
        design=NULL
      } else if (input$find_markers_by == "type"){
        retain <- row.names(subset(metadata, (passage != "normal"  | is.na(passage)) & (source != "cell_line"  | is.na(source)) & (source != "stromal_cells"  | is.na(source))))
        counts_vst <- assay(dds[,retain], "vst")
        selected_meta <- metadata[retain,]
        selected_meta$type <- factor(selected_meta$type)
        selected_meta$source <- factor(selected_meta$source)
        clusters <- factor(selected_meta$type)
        block <- factor(selected_meta$source)
        design <- model.matrix(~block)
      }
      incProgress(1/4)

      pval.type = "any"
      direction =  input$direction
      incProgress(1/4)

      out = findMarkers3(x = counts_vst, clusters = clusters, pval.type="all", direction=direction)
      incProgress(1/4)

      return(list(data=out, retain=retain, selected_meta=selected_meta))
    })
  }) %>%
    bindCache(input$find_markers_by, input$direction)


  makeHeatPlot <- function(){
    retain = markers_data()$retain
    markers_data = markers_data()$data
    selected_meta = markers_data()$selected_meta

    toMatch <- seq(1:input$top_genes)
    markers <- unlist(lapply(1:length(markers_data), function(x) subset(markers_data[[x]], Top %in% toMatch)$Gene))
    markers <- unique(markers)

    output$gene_list <- renderText({ markers })

    if (length(input$column_annotation) > 0){
      expgroups <- as.data.frame(colData(dds)[,input$column_annotation])
      rownames(expgroups) <- colnames(dds)
      colnames(expgroups) <- input$column_annotation
    } else {
      expgroups = NULL
    }

    mydds <- assay(dds[,retain], "vst")

    if(input$remove_source_effect == TRUE){
      mydds <- limma::removeBatchEffect(mydds, batch=selected_meta[["source"]])
    }

    anno_id <- rownames(mydds)
    selectedGenes <- rownames(mydds)[match(markers,rownames(mydds))]
    df <- mydds[selectedGenes, ]

    pheatmap::pheatmap(df,
                       cluster_rows = input$cluster_row,
                       cluster_cols = input$cluster_col,
                       cutree_cols = input$cutree_cols,
                       main=input$title,
                       show_rownames = input$row_labels_hm,
                       show_colnames = input$sample_labels_hm,
                       annotation_col = expgroups,
                       scale=input$scale_by,
                       color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100),
                       fontsize_row = input$fontsize_row)

  }

  output$pheatmap_plot <- renderCachedPlot({

	 makeHeatPlot()
	}, cacheKeyExpr = { list(markers_data(),
			   		input$direction,
			   		input$find_markers_by,
			   		input$top_genes,
					   input$fontsize_row,
					   input$sample_labels_hm,
					   input$row_labels_hm,
					   input$title,
					   input$cutree_cols,
					   input$cluster_row,
					   input$cluster_col,
					   input$remove_source_effect,
					   input$column_annotation,
					   input$color_palette_seq,
					   input$scale_by)})

  callModule(ggDownload, "ggheat", reactive({makeHeatPlot()}))


}




