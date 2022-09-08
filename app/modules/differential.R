
diffInput <- function(id) {

  ns <- NS(id)

  title <- textInput(ns("title"),
                     label = "Plot title",
                     value = "")

  gene_symbol <- selectizeInput(ns("gene_symbol"),
                                "Select gene to plot:",
                                choices = NULL,
                                selected=NULL,
                                multiple=FALSE)


  paste_from_clipboard <- textAreaInput(ns('paste_from_clipboard'), label = 'Include below genes: ', value = "")

  label_de_genes <- checkboxInput(ns('label_de_genes'),"Include DEGs in label?",value = TRUE)

  split_by_type <- checkboxInput(ns('split_by_type'),"Perform differential analysis separately for PDX and primary samples?",value = FALSE)
  keep_intersection <- checkboxInput(ns('keep_intersection'),"Only show genes DE in PDX and primary samples?",value = FALSE)

  multi_gene_symbol <- selectizeInput(ns("multi_gene_symbol"),
                                      "Select genes to plot:",
                                      choices = NULL,
                                      selected=NULL,
                                      multiple=TRUE)



  exprs_values <- selectizeInput(ns("exprs_values"),
                                 "Select expression values:",
                                 choices = c("VST"),
                                 selected="VST",
                                 multiple=FALSE)

  volc_y_axis <- selectizeInput(ns('volc_y_axis'), label = 'Y-axis: ', choices = c("adj.P.Val", "P.Value"), selected="adj.P.Val")

  basemean_cutoff <- numericInput(ns('basemean_cutoff'), label = 'Base mean cutoff: ',value = 5, min = 0,max = 10, step=1)


  top_n_genes <- numericInput(ns('top_n_genes'), label = 'Include top n most variable genes: ', value = 0, min = 0,max = 5000, step=500)

  top_n_from_de_genes <- numericInput(ns('top_n_from_de_genes'), label = 'Show top n most DE genes: ', value = 200, min = 0,max = 8000, step=50)

  sorted_by <- selectizeInput(ns('sorted_by'), label = 'Sorted by: ', choices = c("logFC", "adj.P.Val"), selected="adj.P.Val")

  comp_a <- tags$div(selectizeInput(inputId=ns("comp_a"),
                                    label="Group a",
                                    choices = NULL,
                                    selected=NULL,
                                    multiple=FALSE))

  comp_b <- selectizeInput(inputId=ns("comp_b"),
                           label="Group b",
                           choices = NULL,
                           selected=NULL,
                           multiple=FALSE)

  padj <- numericInput(ns('padj'), label = 'Cutoff for adj. p-value: ', value = 0.001, min = 0.001,max = 1, step=0.01)

  lfc <- sliderInput(ns('lfc'), label = 'Min absolute log2-fold-change required: ', value = 1, min = 0,max = 10, step=1)

  dist_method <- selectInput(inputId=ns("dist_method"),
                             label="Select distance measure",
                             choices = c("euclidean",  "maximum", "manhattan", "canberra", "minkowski"),
                             selected = "euclidean")

  plot_style <- selectInput(ns('plot_style'),
                            "Plot style for gene count",
                            choices = list("dotplot", "boxplot","violin plot"))

  select_all <- actionLink(ns("select_all"), "Select all", style = "color:#18cc7b")
  deselect_all <- actionLink(ns("deselect_all"), "Deselect all", style = "color:#cc3c18")

  group_by <- selectizeInput(inputId=ns("group_by"),
                             label="Group by",
                             choices = NULL,
                             selected=NULL,
                             multiple=TRUE)

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

  deg_res =  textOutput(ns("selected_var"))

  remove_samples <-selectInput(ns('remove_samples'), label = "Remove sample by name:", choices = NULL, selected=NULL, multiple = TRUE)

  include_samples <-selectInput(ns('include_samples'), label = "Include sample by name:", choices = NULL, selected=NULL, multiple = TRUE)


  # list of UIs to be exported to the add
  tagList(title,
#          select_contrasts_for_venn,
       #   exprs_values,
          gene_symbol,
       #   paste_from_clipboard,
          multi_gene_symbol,
#          top_n_genes,
          top_n_from_de_genes,
          volc_y_axis,
		  basemean_cutoff,
          label_de_genes,
          sorted_by,
split_by_type,
keep_intersection,
          fluidRow(column(6,
                          comp_a),
                   column(6, comp_b)
          ),
          deg_res,
          padj,
          lfc,
        #  dist_method,
          group_by,
          plot_style,
          facet_by_col,
          facet_by_row,
          hr(),
          remove_batch_effect,
         hr(),
helpText("Differential expression analysis is performed using limma on voom normalized counts.")
		  )
}


diffOutput <- function(id) {
  ## Define the different types of plots that will be shown
  ns <- NS(id)


  tagList(
    ## cluster
    tabsetPanel(id = ns("tabs"),
                tabPanel("Heatmap",
                         shinydashboard::box( width = 12, title = "Gene expression heatmap",
                                              fluidRow(column(width = 3,
                                                              checkboxInput(ns('sample_labels_hm'),
                                                                            "Display sample labels?",
                                                                            value = TRUE)),
                                                       column(width = 3,
                                                              checkboxInput(ns('row_labels_hm'),
                                                                            "Display gene IDs?",
                                                                            value = FALSE)),
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
                                                                          choices = c(row.names(subset(RColorBrewer:::brewer.pal.info, category=="div"  | category=="seq")),c('viridis', 'inferno', 'plasma', 'magma')),
                                                                          selected = "RdYlBu" ,multiple = FALSE))
                                              ),

                                              fluidRow(
                                                column(width = 3,
                                                       numericInput(ns('fontsize_row'), label = 'Row font size:', value = 10, min = 1,max = 20, step=1 )),
                                                column(width = 3,
                                                       numericInput(ns('cell_width'), label = 'Cell width in points. If 0, then auto:',  value = 0, min = 0,max = 100, step=2)),
                                                column(width = 3,
                                                       numericInput(ns('cell_height'), label = 'Cell height in points. If 0, then auto:',value = 0, min = 0,max = 20, step=0.10)),
                                                column(width = 3,
                                                       numericInput(ns('breaks'), label = 'Max scale range. If 0, then auto:',value = 0, min = 0,max = 10, step=0.50))
                                              ),
                                             checkboxInput(ns('keep_row_order'),"Keep row order the same?",value = FALSE),

                                             actionLink(ns("toggle_venn_scatter"), "Hide/show venn diagram / logFC scatter"),
                                             br(),
                                              downloadLink(ns("downloadData"), "Download gene table"),
                                             plotOutput(ns("venn_plot"), height="300px"),
                                              plotOutput(ns("pheatmap_plot"), height="900px"),
                                             ggDownloadUI(ns("ggheat")),
                                              ),

                         collapsible = FALSE
                ),
                tabPanel("Volcano plot",
                         shinydashboard::box( width = 12, title = "Volcano plot",
                                              fluidRow(column(width = 2,
                                                              numericInput(ns('pointSize'),
                                                                           "Point size:",value = 2, min = 1,max = 20, step=1)),
                                                       column(width = 2,
                                                              numericInput(ns('labSize'),
                                                                           "Label size:", value = 3, min = 1,max = 20, step=1)),
                                                       column(width = 2,
                                                              numericInput(ns('colAlpha'),
                                                                           "Point transparency:", value = 0.5, min = 0.1,max = 1, step=0.1)),
                                                       column(width=2,
                                                              numericInput(ns('shape'),
                                                                           "Point shape:",value = 19, min = 0,max = 25, step=1)),
                                                       column(width=4,
                                                              selectInput(ns('col'),
                                                                          "Point color (4):",choices = colors(), selected =c('goldenrod1', 'black', 'tomato2', 'maroon4'),
                                                                          multiple = TRUE)),
                                                       # column(width=3,
                                                       #        selectInput(ns('shape'),
                                                       #                    "Point Shape:",choices = seq(0,25), selected = 19,
                                                       #                    multiple = TRUE))
                                              ),
                                              fluidRow(column(width=3,
                                                              selectInput(ns('labFace'),
                                                                          "Label face:",choices = c('plain', 'bold', 'italic', 'bold-italic'), selected = "plain",
                                                                          multiple = FALSE)),
                                                       column(width=3,
                                                              selectInput(ns('legendPosition'),
                                                                          "Legend position:",choices = c('top', 'bottom', 'left', 'right'), selected = "bottom",
                                                                          multiple = FALSE)),
                                                       column(width = 3,
                                                              numericInput(ns('legendIconSize'),
                                                                           "Legend icon size:", value = 4, min = 1,max = 20, step=1)),
                                                       column(width = 3,
                                                              numericInput(ns('legendLabSize'),
                                                                           "Legend label size:", value = 14, min = 1,max = 20, step=1))

                                              ),
                                              fluidRow(column(width = 3,
                                                              checkboxInput(ns('boxedLabels'),
                                                                            "Boxed labels?",
                                                                            value = FALSE)),
                                                       column(width = 3,
                                                              checkboxInput(ns('drawConnectors'),
                                                                            "Draw connectors?",
                                                                            value = FALSE))

                                              ),
                                              plotOutput(ns("volcano_plot"), height="850px"),
											helpText("To change color of points, you must select four colors corresponding to < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff."),
											ggDownloadUI(ns("ggvol")),

                                              collapsible = FALSE)
                ),

                tabPanel("MA plot",
                         shinydashboard::box( width = 12, title = "MA plot",
                                              fluidRow(column(width = 2,
                                                              numericInput(ns('ma_pointSize'),
                                                                           "Point size:",value = 2, min = 1,max = 20, step=1)),
                                                       column(width = 2,
                                                              numericInput(ns('ma_labSize'),
                                                                           "Label size:", value = 3, min = 1,max = 20, step=1)),
                                                       column(width = 2,
                                                              numericInput(ns('ma_colAlpha'),
                                                                           "Point transparency:", value = 0.8, min = 0.1,max = 1, step=0.1)),
                                                       column(width=2,
                                                              numericInput(ns('ma_shape'),
                                                                           "Point shape:",value = 19, min = 0,max = 25, step=1)),
                                                       column(width=4,
                                                              selectInput(ns('ma_col'),
                                                              "Point color (4):",choices = colors(), selected =c('goldenrod1', 'black', 'tomato2', 'maroon4'),
                                                                          multiple = TRUE)),
                                              ),
                                              fluidRow(column(width=3,
                                                              selectInput(ns('ma_labFace'),
                                                                          "Label face:",choices = c('plain', 'bold', 'italic', 'bold-italic'), selected = "plain",
                                                                          multiple = FALSE)),
                                                       column(width=3,
                                                              selectInput(ns('ma_legendPosition'),
                                                                          "Legend position:",choices = c('top', 'bottom', 'left', 'right'), selected = "bottom",
                                                                          multiple = FALSE)),
                                                       column(width = 3,
                                                              numericInput(ns('ma_legendIconSize'),
                                                                           "Legend icon size:", value = 4, min = 1,max = 20, step=1)),
                                                       column(width = 3,
                                                              numericInput(ns('ma_legendLabSize'),
                                                                           "Legend label size:", value = 14, min = 1,max = 20, step=1))

                                              ),
                                              fluidRow(column(width = 3,
                                                              checkboxInput(ns('ma_boxedLabels'),
                                                                            "Boxed labels?",
                                                                            value = FALSE)),
                                                       column(width = 3,
                                                              checkboxInput(ns('ma_drawConnectors'),
                                                                            "Draw connectors?",
                                                                            value = FALSE))

                                              ),
                                              plotOutput(ns("ma_plot"), height="850px"),
                                              ggDownloadUI(ns("ggma")),
                                              collapsible = FALSE)
                ),

                tabPanel("Tables",
                         shinydashboard::box( width = 12, title = "Tables",
			                                      DT::dataTableOutput(ns('deg_table'), height="1000px"),
                                              collapsible = FALSE)
                )

    )

  )

}


diffMod <- function(input, output, session, dds = dds, metadata, fit = fit, fit_pdx=fit_pdx, fit_primary=fit_primary) {
 # shinyjs::toggle("venn_plot")
	updateSelectizeInput(session, "multi_gene_symbol", choices = c(row.names(dds)), selected=NULL, server=TRUE)
	updateSelectInput(session, "group_by", choices = c("type"), selected = "type")
	updateSelectInput(session, "column_annotation", choices = c("", names(metadata)), selected="type")
	updateSelectInput(session, "remove_batch_effect", choices = c("none","source", "batch", "patient", "purity"), selected=NULL)
	updateSelectInput(session, "comp_a", choices = make.names(levels(dds$type_design)) ,selected=make.names(levels(dds$type_design))[1])
	updateSelectInput(session, "comp_b", choices = make.names(levels(dds$type_design)) ,selected=make.names(levels(dds$type_design))[2])
	updateSelectInput(session, "col", choices = colors(), selected=c('goldenrod1', 'black', 'tomato2', 'maroon4'))


	observeEvent(input$split_by_type, {
	  if(input$split_by_type == TRUE){
	    common_sets = intersect(colnames(fit_pdx$coefficients), colnames(fit_primary$coefficients))
	    updateSelectInput(session, "comp_a", choices = make.names(common_sets) ,selected=make.names(common_sets)[1])
	    updateSelectInput(session, "comp_b", choices = make.names(common_sets) ,selected=make.names(common_sets)[2])
	  } else {
	    updateSelectInput(session, "comp_a", choices = make.names(levels(dds$type_design)) ,selected=make.names(levels(dds$type_design))[1])
	    updateSelectInput(session, "comp_b", choices = make.names(levels(dds$type_design)) ,selected=make.names(levels(dds$type_design))[2])
	  }
	})


	observeEvent(input$split_by_type, {
	  if(input$split_by_type == TRUE){
	    shinyjs::show("keep_intersection")
	    shinyjs::show("toggle_venn_scatter")
	    shinyjs::show("keep_row_order")
	  } else {
	    shinyjs::hide("keep_intersection")
	    shinyjs::hide("toggle_venn_scatter")
	    shinyjs::hide("keep_row_order")
	    shinyjs::hide("venn_plot")
	  }

	})

  observeEvent(input$tabs, {
    if (input$tabs ==  "Heatmap"){
      shinyjs::show("top_n_from_de_genes")
      shinyjs::show("comp_b")
      shinyjs::show("comp_a")
      shinyjs::show("lfc")
      shinyjs::show("padj")
      shinyjs::show("remove_batch_effect")
      shinyjs::hide("multi_gene_symbol")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("group_by")
      shinyjs::hide("plot_style")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("show")
     # shinyjs::hide("paste_from_clipboard")
      shinyjs::show("sorted_by")
      shinyjs::show("deg_res")
      shinyjs::show("select_samples_by")
      shinyjs::show("selected_samples")
      shinyjs::show("select_all")
      shinyjs::show("deselect_all")
      shinyjs::show("remove_samples")
      shinyjs::show("include_samples")
      shinyjs::hide("label_de_genes")
      shinyjs::hide("volc_y_axis")
	  shinyjs::hide("basemean_cutoff")
    } else if (input$tabs ==  "Volcano plot"){
      shinyjs::hide("top_n_genes")
      shinyjs::hide("top_n_from_de_genes")
      shinyjs::show("label_de_genes")
      shinyjs::show("comp_b")
      shinyjs::show("comp_a")
      shinyjs::show("lfc")
      shinyjs::show("padj")
      shinyjs::hide("remove_batch_effect")
      shinyjs::show("multi_gene_symbol")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("group_by")
      shinyjs::hide("plot_style")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("show")
	  shinyjs::hide("title")
    #  shinyjs::show("paste_from_clipboard")
      shinyjs::hide("sorted_by")
      shinyjs::show("deg_res")
      shinyjs::hide("select_samples_by")
      shinyjs::hide("selected_samples")
      shinyjs::hide("select_all")
      shinyjs::hide("deselect_all")
      shinyjs::hide("remove_samples")
      shinyjs::hide("include_samples")
      shinyjs::show("volc_y_axis")
	  shinyjs::hide("basemean_cutoff")
    } else if (input$tabs ==  "MA plot"){
      shinyjs::hide("top_n_genes")
      shinyjs::hide("top_n_from_de_genes")
      shinyjs::show("label_de_genes")
      shinyjs::show("comp_b")
      shinyjs::show("comp_a")
      shinyjs::show("lfc")
      shinyjs::show("padj")
      shinyjs::hide("remove_batch_effect")
      shinyjs::show("multi_gene_symbol")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("group_by")
      shinyjs::hide("plot_style")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("show")
	  shinyjs::hide("title")
      shinyjs::hide("sorted_by")
      shinyjs::show("deg_res")
      shinyjs::hide("select_samples_by")
      shinyjs::hide("selected_samples")
      shinyjs::hide("select_all")
      shinyjs::hide("deselect_all")
      shinyjs::hide("remove_samples")
      shinyjs::hide("include_samples")
      shinyjs::hide("volc_y_axis")
	  shinyjs::show("basemean_cutoff")
    } else if (input$tabs ==  "Tables"){
      shinyjs::hide("top_n_genes")
      shinyjs::hide("top_n_from_de_genes")
      shinyjs::hide("label_de_genes")
      shinyjs::show("comp_b")
      shinyjs::show("comp_a")
      shinyjs::show("lfc")
      shinyjs::show("padj")
      shinyjs::hide("remove_batch_effect")
      shinyjs::hide("multi_gene_symbol")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("group_by")
      shinyjs::hide("plot_style")
      shinyjs::hide("facet_by_row")
      shinyjs::hide("facet_by_col")
      shinyjs::hide("show")
	  shinyjs::hide("title")
      shinyjs::hide("sorted_by")
      shinyjs::show("deg_res")
      shinyjs::hide("select_samples_by")
      shinyjs::hide("selected_samples")
      shinyjs::hide("select_all")
      shinyjs::hide("deselect_all")
      shinyjs::hide("remove_samples")
      shinyjs::hide("include_samples")
      shinyjs::hide("volc_y_axis")
	  shinyjs::hide("basemean_cutoff")
    }
  })

  timer <- reactive(list(comp_a = input$comp_a,
						 comp_b = input$comp_b,
						 padj = input$padj,
						 lfc = input$lfc,
                         top_n_from_de_genes=input$top_n_from_de_genes)) %>% debounce(0)



  observeEvent(input$toggle_venn_scatter,{
    shinyjs::toggle("venn_plot")
  })


  degs <- reactive({
    shiny::validate(
      need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
           message = "Select two different groups to compare against.")
    )
    withProgress(message = 'Performing differential expression analysis...', {
      incProgress(0)

      if(input$split_by_type == FALSE){
        myargs = list(paste(input$comp_a,input$comp_b,sep="-"), levels=fit$design)
        contrast.matrix = do.call(makeContrasts, myargs)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        degs <- topTable(fit2, coef=1, n=Inf, adjust="BH")
      } else {

        shiny::validate(
          need(make.names(timer()$comp_a) %in% colnames(fit_pdx$design),
               message = "")
        )
        shiny::validate(
          need(make.names(timer()$comp_b) %in% colnames(fit_pdx$design),
               message = "")
        )
        shiny::validate(
          need(make.names(timer()$comp_a) %in% colnames(fit_primary$design),
               message = "")
        )
        shiny::validate(
          need(make.names(timer()$comp_b) %in% colnames(fit_primary$design),
               message = "")
        )



        myargs = list(paste(input$comp_a,input$comp_b,sep="-"), levels=fit_pdx$design)
        contrast.matrix = do.call(makeContrasts, myargs)
        fit2_pdx <- contrasts.fit(fit_pdx, contrast.matrix)
        fit2_pdx <- eBayes(fit2_pdx)
        degs_pdx <- topTable(fit2_pdx, coef=1, n=Inf, adjust="BH")

        myargs = list(paste(input$comp_a,input$comp_b,sep="-"), levels=fit_primary$design)
        contrast.matrix = do.call(makeContrasts, myargs)
        fit2_primary <- contrasts.fit(fit_primary, contrast.matrix)
        fit2_primary <- eBayes(fit2_primary)
        degs_primary <- topTable(fit2_primary, coef=1, n=Inf, adjust="BH")
        degs <- list(pdx=degs_pdx, primary=degs_primary)
      }
      return(degs)
    })
  }) %>%
    bindCache(timer()$comp_a, timer()$comp_b, timer()$padj, input$split_by_type)


    output$selected_var <- renderText({
      if(input$split_by_type == FALSE){
        paste(nrow(subset(degs(), adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc)), "genes detected as DE: ", paste(timer()$comp_a," - ",timer()$comp_b), paste("adj. p < ", timer()$padj, "; abs(LFC) > ", timer()$lfc))
      }
      })


    output$venn_plot <- renderCachedPlot({
      shiny::validate(
        need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
             message = "Select two different groups to compare against.")
      )

      if(input$split_by_type == TRUE){
        primary_degs = degs()$primary
        primary_degs <-  subset(primary_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()

        pdx_degs = degs()$pdx
        pdx_degs <-  subset(pdx_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()

        ll <- list(PDX=row.names(pdx_degs),primary=row.names(primary_degs))

        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

        venn.plot <- VennDiagram::venn.diagram(ll,
                                               alpha = c(0.3, 0.3),
                                               fill=c("dodgerblue", "goldenrod1"),
                                               cex = 1,cat.fontface = 2, filename = NULL, margin = 0.1)#, fill = c("dodgerblue", "goldenrod1", "red", "green"))

        grid::grid.draw(venn.plot);

        primary_degs_unfilt <- degs()$primary %>% as.data.frame()
        primary_degs_unfilt$sig <- ifelse(primary_degs_unfilt$adj.P.Val < input$padj & abs(primary_degs_unfilt$logFC) > input$lfc, TRUE, FALSE)


        pdx_degs_unfilt <- degs()$pdx %>% as.data.frame()
        pdx_degs_unfilt$sig <- ifelse(pdx_degs_unfilt$adj.P.Val < input$padj & abs(pdx_degs_unfilt$logFC) > input$lfc, TRUE, FALSE)

        merged_dfs = merge(primary_degs_unfilt, pdx_degs_unfilt, by="row.names")
        merged_dfs$color <- "grey"
        merged_dfs$color <- ifelse(merged_dfs$logFC.x > 0 & merged_dfs$logFC.y > 0, "red", ifelse(merged_dfs$logFC.x < 0 & merged_dfs$logFC.y < 0, "blue", "grey"))


      p2 <-   ggplot(merged_dfs, aes(x=logFC.x, y=logFC.y, color=color, label=Row.names)) + geom_point(size=1, alpha=0.4) + theme_bw() +
          scale_color_manual(values=c("blue", "grey", "red"),  labels=c( "different", "negative", "positive")) +
          xlim(floor( min(c(merged_dfs$logFC.x, merged_dfs$logFC.y))), ceiling( max(c(merged_dfs$logFC.x, merged_dfs$logFC.y))))+
          ylim(floor( min(c(merged_dfs$logFC.x, merged_dfs$logFC.y))), ceiling( max(c(merged_dfs$logFC.x, merged_dfs$logFC.y)))) +
          theme(legend.position="bottom") +  guides(color=guide_legend(title="fold-change:")) +   geom_smooth(method=lm, se=F, color="black") +
          xlab(paste0("primary logFC : ", timer()$comp_a, " - ", timer()$comp_b)) +
          ylab(paste0("PDX logFC : ", timer()$comp_a, " - ", timer()$comp_b))

      gridExtra::grid.arrange(grobs = list(grid::gTree(children=venn.plot),
                                           p2),
                              ncol = 2, widths=c(1,3))


      }}, cacheKeyExpr = { list(degs(),
                             timer()$comp_a,
                             timer()$comp_b,
                             timer()$padj,
                             timer()$lfc,
                             timer()$top_n_from_de_genes,
                             input$split_by_type)})


    makeHeatPlot <- function(){

      shiny::validate(
        need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
             message = "Select two different groups to compare against.")
      )


      if(input$split_by_type == FALSE){
        coldata = data.frame(colData(dds))
        coldata$type_design <- make.names(coldata$type_design)
        samps = subset(coldata, type_design == make.names(timer()$comp_a) | type_design == make.names(timer()$comp_b))

        shiny::validate(
          need(timer()$top_n_from_de_genes > 0,
               message = "Number of genes to plot must be greater than 0.")
        )

        object = dds[,row.names(samps)]

        vst_data <- assay(object, "vst")

        if (input$remove_batch_effect != "none"){
          tryCatch(vst_data <- limma::removeBatchEffect(vst_data, batch=colData(object)[[input$remove_batch_effect]]), error=function(e) vst_data <- assay(object, "vst"))
        }


        degs = degs()
        degs <-  subset(degs , adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        if(input$sorted_by == "adj.P.Val"){
          degs <- head(degs[order(degs$adj.P.Val),], timer()$top_n_from_de_genes)
        } else {
          degs <- head(degs[order(abs(degs$logFC), decreasing=T),], timer()$top_n_from_de_genes)
        }
        degs <- row.names(degs)

        df <- vst_data[degs, ]

        output$downloadData <- downloadHandler(
          filename = function() {
            paste("hm_table", ".csv", sep = "")
          },
          content = function(file) {
            write.csv(sample_data(), file, row.names = TRUE)
          })

        sample_data <- reactive({as.data.frame(df)})

        if (length(input$column_annotation) > 0){
          expgroups <- as.data.frame(colData(dds)[,input$column_annotation])
          rownames(expgroups) <- colnames(dds)
          colnames(expgroups) <- input$column_annotation
        } else {
          expgroups = NULL
        }


        if(input$cell_height == 0) {cellheight = NA}  else {cellheight = input$cell_height }
        if(input$cell_width == 0) {cellwidth = NA}  else { cellwidth = input$cell_width }

        shiny::validate(
          need(length(degs) > 2,
               message = paste("Not enough genes after filtering to draw a heatmap.", "Only the following genes selected:", paste(degs, collapse=", ")))
        )


        if(input$scale_by !="none"){
          df <- df[apply(df, MARGIN = 1, FUN = function(x) sd(x) != 0),]
          shiny::validate(
            need( nrow(df) > 2,
                  message = paste("Turn off row scaling to view selected genes.", "The following genes were selected:", paste(degs, collapse=", ")))
          )
        }

        if(input$color_palette_seq %in% c('viridis', 'inferno', 'plasma', 'magma')){
          x=input$color_palette_seq
          require('viridis')
          viridis_ramp <- get(x, envir = asNamespace("viridis"), mode = "function")
          np <- 100
          hm_color =  viridis_ramp(np)
        } else {
          hm_color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100)
        }

        if(input$breaks == 0){
          breaks = NA
        } else {
          breaks = seq(-(input$breaks), input$breaks, length.out = 100)
        }

        ph <-  pheatmap::pheatmap(df,  cluster_rows = input$cluster_row,
                                  cluster_cols = input$cluster_col,
                                  cutree_cols = input$cutree_cols,
                                  main=input$title,
                                  fontsize_row = input$fontsize_row,
                                  show_rownames = input$row_labels_hm ,
                                  show_colnames = input$sample_labels_hm,
                                  annotation_col = expgroups,
                                  scale=input$scale_by,
                                  cellheight = cellheight, cellwidth = cellwidth,
                                  color = hm_color, breaks=breaks)
      } else {


        if(input$breaks == 0){
          breaks = NA
        } else {
          breaks = seq(-(input$breaks), input$breaks, length.out = 100)
        }


        output$downloadData <- downloadHandler(
          filename = function() {
            paste("hm_table", ".csv", sep = "")
          },
          content = function(file) {
            write.csv(sample_data(), file, row.names = TRUE)
          })


        if(input$keep_intersection == TRUE){
          primary_degs = degs()$primary
          primary_degs <-  subset(primary_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
          # if(input$sorted_by == "adj.P.Val"){
          #   primary_degs <- head(primary_degs[order(primary_degs$adj.P.Val),], input$top_n_from_de_genes)
          # } else {
          #   primary_degs <- head(primary_degs[order(abs(primary_degs$logFC), decreasing=T),], input$top_n_from_de_genes)
          # }
          primary_degs <- row.names(primary_degs)

          pdx_degs = degs()$pdx
          pdx_degs <-  subset(pdx_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
          # if(input$sorted_by == "adj.P.Val"){
          #   pdx_degs <- head(pdx_degs[order(pdx_degs$adj.P.Val),], input$top_n_from_de_genes)
          # } else {
          #   pdx_degs <- head(pdx_degs[order(abs(pdx_degs$logFC), decreasing=T),], input$top_n_from_de_genes)
          # }
          pdx_degs <- row.names(pdx_degs)

          inter <- intersect(primary_degs, pdx_degs)
          inter <- head(inter, input$top_n_from_de_genes)
          shiny::validate(
            need(length(inter) > 2,
                 message = paste("Not enough genes after filtering to draw a heatmap.", "Only the following genes in intersection:", paste(inter, collapse=", ")))
          )
        }



        coldata = data.frame(colData(dds))
        coldata$type_design <- make.names(coldata$type_design)
        coldata = subset(coldata, source == "PDX")
        samps = subset(coldata, type_design == make.names(input$comp_a) | type_design == make.names(input$comp_b))
        object = dds[,row.names(samps)]
        vst_data <- assay(object, "vst")

        pdx_degs = degs()$pdx
        pdx_degs <-  subset(pdx_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        if(input$sorted_by == "adj.P.Val"){
          pdx_degs <- head(pdx_degs[order(pdx_degs$adj.P.Val),], input$top_n_from_de_genes)
        } else {
          pdx_degs <- head(pdx_degs[order(abs(pdx_degs$logFC), decreasing=T),], input$top_n_from_de_genes)
        }
        pdx_degs <- row.names(pdx_degs)

        if(input$keep_intersection == TRUE){ pdx_degs = inter }

        df_pdx <- vst_data[pdx_degs, ]


        if (length(input$column_annotation) > 0){
          expgroups <- as.data.frame(samps[,input$column_annotation])
          rownames(expgroups) <- colnames(vst_data)
          colnames(expgroups) <- input$column_annotation
        } else {
          expgroups = NULL
        }


        if(input$cell_height == 0) {cellheight = NA}  else {cellheight = input$cell_height }
        if(input$cell_width == 0) {cellwidth = NA}  else { cellwidth = input$cell_width }

        shiny::validate(
          need(length(pdx_degs) > 2,
               message = paste("Not enough genes after filtering to draw a heatmap.", "Only the following genes selected:", paste(pdx_degs, collapse=", ")))
        )



        if(input$scale_by !="none"){
          df_pdx <- df_pdx[apply(df_pdx, MARGIN = 1, FUN = function(x) sd(x) != 0),]
          shiny::validate(
            need( nrow(df_pdx) > 2,
                  message = paste("Turn off row scaling to view selected genes.", "The following genes were selected:", paste(pdx_degs, collapse=", ")))
          )
        }

        if(input$color_palette_seq %in% c('viridis', 'inferno', 'plasma', 'magma')){
          x=input$color_palette_seq
          require('viridis')
          viridis_ramp <- get(x, envir = asNamespace("viridis"), mode = "function")
          np <- 100
          hm_color =  viridis_ramp(np)
        } else {
          hm_color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100)
        }

        ph <-  pheatmap::pheatmap(df_pdx,  cluster_rows = input$cluster_row,
                                  cluster_cols = input$cluster_col,
                                  cutree_cols = input$cutree_cols,
                                  main="PDX",
                                  fontsize_row = input$fontsize_row,
                                  show_rownames = input$row_labels_hm ,
                                  show_colnames = input$sample_labels_hm,
                                  annotation_col = expgroups,
                                  scale=input$scale_by,
                                  cellheight = cellheight, cellwidth = cellwidth,
                                  color = hm_color, silent=T, breaks=breaks)



        # primary
        coldata = data.frame(colData(dds))
        coldata$type_design <- make.names(coldata$type_design)
        coldata = subset(coldata, source == "primary")
        samps = subset(coldata, type_design == make.names(input$comp_a) | type_design == make.names(input$comp_b))
        object = dds[,row.names(samps)]
        vst_data <- assay(object, "vst")

        primary_degs = degs()$primary
        primary_degs <-  subset(primary_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        if(input$sorted_by == "adj.P.Val"){
          primary_degs <- head(primary_degs[order(primary_degs$adj.P.Val),], input$top_n_from_de_genes)
        } else {
          primary_degs <- head(primary_degs[order(abs(primary_degs$logFC), decreasing=T),], input$top_n_from_de_genes)
        }
        primary_degs <- row.names(primary_degs)
        if(input$keep_intersection == TRUE){ primary_degs = inter }
        df_primary <- vst_data[primary_degs, ]


        if (length(input$column_annotation) > 0){
          expgroups <- as.data.frame(samps[,input$column_annotation])
          rownames(expgroups) <- colnames(vst_data)
          colnames(expgroups) <- input$column_annotation
        } else {
          expgroups = NULL
        }


        if(input$cell_height == 0) {cellheight = NA}  else {cellheight = input$cell_height }
        if(input$cell_width == 0) {cellwidth = NA}  else { cellwidth = input$cell_width }

        shiny::validate(
          need(length(primary_degs) > 2,
               message = paste("Not enough genes after filtering to draw a heatmap.", "Only the following genes selected:", paste(primary_degs, collapse=", ")))
        )


        if(input$scale_by !="none"){
          df_primary <- df_primary[apply(df_primary, MARGIN = 1, FUN = function(x) sd(x) != 0),]
          shiny::validate(
            need( nrow(df_primary) > 2,
                  message = paste("Turn off row scaling to view selected genes.", "The following genes were selected:", paste(primary_degs, collapse=", ")))
          )
        }

        if(input$color_palette_seq %in% c('viridis', 'inferno', 'plasma', 'magma')){
          x=input$color_palette_seq
          require('viridis')
          viridis_ramp <- get(x, envir = asNamespace("viridis"), mode = "function")
          np <- 100
          hm_color =  viridis_ramp(np)
        } else {
          hm_color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100)
        }

        joined = merge(as.data.frame(df_pdx), as.data.frame(df_primary), by="row.names")
        sample_data <- reactive({as.data.frame(joined)})


        row.order = ph$tree_row$order
        row.order = row.names(df_pdx[row.order,])


        if(input$keep_intersection == TRUE & input$keep_row_order == TRUE){

          ph2 <-  pheatmap::pheatmap(df_primary[row.order,],  cluster_rows = F,
                                     cluster_cols = input$cluster_col,
                                     cutree_cols = input$cutree_cols,
                                     main="Primary",
                                     fontsize_row = input$fontsize_row,
                                     show_rownames = input$row_labels_hm ,
                                     show_colnames = input$sample_labels_hm,
                                     annotation_col = expgroups,
                                     scale=input$scale_by,
                                     cellheight = cellheight, cellwidth = cellwidth,
                                     color = hm_color, silent=T, breaks=breaks)
        } else {

          ph2 <-  pheatmap::pheatmap(df_primary,  cluster_rows = input$cluster_row,
                                     cluster_cols = input$cluster_col,
                                     cutree_cols = input$cutree_cols,
                                     main="Primary",
                                     fontsize_row = input$fontsize_row,
                                     show_rownames = input$row_labels_hm ,
                                     show_colnames = input$sample_labels_hm,
                                     annotation_col = expgroups,
                                     scale=input$scale_by,
                                     cellheight = cellheight, cellwidth = cellwidth,
                                     color = hm_color, silent=T, breaks=breaks)
        }


        gridExtra::grid.arrange(grobs = list(ph[[4]],
                                             ph2[[4]]),
                                ncol = 2)
      }


      }

  output$pheatmap_plot <- renderCachedPlot({
      shiny::validate(
        need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
             message = "Select two different groups to compare against.")
      )


    makeHeatPlot()

  }, cacheKeyExpr = { list(degs(),
                           timer()$comp_a,
                           timer()$comp_b,
                           timer()$padj,
                           timer()$lfc,
                           timer()$top_n_from_de_genes,
                           input$sorted_by,
						                input$scale_by,
						                input$sample_labels_hm,
						                input$row_labels_hm,
						                input$fontsize_row,
						                input$title,
						                input$cutree_cols,
						                input$cluster_col,
						                input$cluster_row,
						                input$color_palette_seq,
						                input$scale_by,
						                input$cell_height,
						                input$cell_width,
						                input$column_annotation,
						                input$remove_batch_effect,
						                input$split_by_type,
						                input$keep_intersection,
						                input$keep_row_order,
						                input$breaks)})



  makeVolPlot <- function() {

    shiny::validate(
      need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
           message = "Select two different groups to compare against.")
    )


    shiny::validate(
      need(length(input$col) == 4,
           message = "Select exactly four colors corresponding to < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff.")
    )


    if(input$split_by_type == FALSE){
      degs_all = degs()

      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }

      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      pv <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                             lab = rownames(degs_all),
                                             x = 'logFC',
                                             y = input$volc_y_axis,
                                             pCutoff = timer()$padj,
                                             FCcutoff = timer()$lfc,
                                             # title=input$title,
                                             selectLab = selectedGenes,
                                             pointSize=input$pointSize,
                                             labSize=input$labSize,
                                             colAlpha=input$colAlpha,
                                             shape=as.numeric(input$shape),
                                             legendPosition=input$legendPosition,
                                             subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                             legendLabSize=input$legendLabSize,
                                             legendIconSize=input$legendIconSize,
                                             labFace=input$labFace,
                                             boxedLabels=input$boxedLabels,
                                             drawConnectors=input$drawConnectors,
                                             gridlines.major = FALSE,
                                             gridlines.minor = FALSE,
                                             widthConnectors = 0.25,
                                             col=input$col)
      pv

    } else {



      if(input$keep_intersection == TRUE){
        primary_degs = degs()$primary
        primary_degs <-  subset(primary_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        primary_degs <- row.names(primary_degs)
        pdx_degs = degs()$pdx
        pdx_degs <-  subset(pdx_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        pdx_degs <- row.names(pdx_degs)
        inter <- intersect(primary_degs, pdx_degs)
      }



      #  pdx
      degs_all = degs()$pdx

      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }

      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      if(input$keep_intersection == TRUE){ selectedGenes = inter }

      pv <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                             lab = rownames(degs_all),
                                             x = 'logFC',
                                             y = input$volc_y_axis,
                                             pCutoff = timer()$padj,
                                             FCcutoff = timer()$lfc,
                                             title="Volcano plot -- PDX samples",
                                             selectLab = selectedGenes,
                                             pointSize=input$pointSize,
                                             labSize=input$labSize,
                                             colAlpha=input$colAlpha,
                                             shape=as.numeric(input$shape),
                                             legendPosition=input$legendPosition,
                                             subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                             legendLabSize=input$legendLabSize,
                                             legendIconSize=input$legendIconSize,
                                             labFace=input$labFace,
                                             boxedLabels=input$boxedLabels,
                                             drawConnectors=input$drawConnectors,
                                             gridlines.major = FALSE,
                                             gridlines.minor = FALSE,
                                             widthConnectors = 0.25,
                                             col=input$col)


      #  primary
      degs_all = degs()$primary

      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }

      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      if(input$keep_intersection == TRUE){ selectedGenes = inter }

      pv2 <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                              lab = rownames(degs_all),
                                              x = 'logFC',
                                              y = input$volc_y_axis,
                                              pCutoff = timer()$padj,
                                              FCcutoff = timer()$lfc,
                                              title="Volcano plot -- primary samples",
                                              selectLab = selectedGenes,
                                              pointSize=input$pointSize,
                                              labSize=input$labSize,
                                              colAlpha=input$colAlpha,
                                              shape=as.numeric(input$shape),
                                              legendPosition=input$legendPosition,
                                              subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                              legendLabSize=input$legendLabSize,
                                              legendIconSize=input$legendIconSize,
                                              labFace=input$labFace,
                                              boxedLabels=input$boxedLabels,
                                              drawConnectors=input$drawConnectors,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              widthConnectors = 0.25,
                                              col=input$col)

      grid.arrange(pv, pv2, nrow = 2)
    }




  }

  output$volcano_plot <- renderCachedPlot({

      shiny::validate(
        need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
             message = "Select two different groups to compare against.")
      )


      shiny::validate(
        need(length(input$col) == 4,
             message = "Select exactly four colors corresponding to < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff.")
			 )

      makeVolPlot()

  }, cacheKeyExpr = { list(degs(),
                           timer()$comp_a,
                           timer()$comp_b,
                           timer()$padj,
                           timer()$lfc,
                           timer()$top_n_from_de_genes,
                           input$col,
                           input$label_de_genes,
                           input$multi_gene_symbol,
                           input$drawConnectors,
                           input$boxedLabels,
                           input$labFace,
                           input$legendIconSize,
                           input$legendLabSize,
                           input$legendPosition,
                           input$shape,
                           input$colAlpha,
                           input$labSize,
                           input$pointSize,
                           input$volc_y_axis,
                           input$split_by_type,
                           input$keep_intersection)})





  makeMaPlot <- function(){

    shiny::validate(
      need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
           message = "Select two different groups to compare against.")
    )


    shiny::validate(
      need(length(input$ma_col) == 4,
           message = "Select exactly four colors corresponding to < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff.")
    )


    if(input$split_by_type == FALSE){
      degs_all = degs()


      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }


      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      degs_all$baseMeanNew <- 1 / (10^log((2^degs_all$AveExpr) + 1))

      pv <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                             lab = rownames(degs_all),
                                             x = 'logFC',
                                             y = "baseMeanNew", #input$volc_y_axis,
                                             pCutoff =  1/10^as.numeric(input$basemean_cutoff),
                                             FCcutoff = timer()$lfc,
                                             xlab = bquote(~Log[2]~ 'fold change'),
                                             ylab = bquote(~Log[e]~ 'base mean + 1'),
                                             legendLabels = c('NS', expression(Log[2]~FC), 'Mean expr.', expression(Mean-expr~and~log[2]~FC)),
                                             title="MA plot",
                                             selectLab = selectedGenes,
                                             pointSize=input$ma_pointSize,
                                             labSize=input$ma_labSize,
                                             colAlpha=input$ma_colAlpha,
                                             shape=as.numeric(input$ma_shape),
                                             legendPosition=input$ma_legendPosition,
                                             subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                             legendLabSize=input$ma_legendLabSize,
                                             legendIconSize=input$ma_legendIconSize,
                                             labFace=input$ma_labFace,
                                             boxedLabels=input$ma_boxedLabels,
                                             drawConnectors=input$ma_drawConnectors,
                                             gridlines.major = FALSE,
                                             gridlines.minor = FALSE,
                                             widthConnectors = 0.25,
                                             ylim = c(0, 10),
                                             col=input$ma_col) + coord_flip()
      pv

    } else {

      if(input$keep_intersection == TRUE){
        primary_degs = degs()$primary
        primary_degs <-  subset(primary_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        primary_degs <- row.names(primary_degs)
        pdx_degs = degs()$pdx
        pdx_degs <-  subset(pdx_degs, adj.P.Val < input$padj  & abs(logFC) > input$lfc) %>% as.data.frame()
        pdx_degs <- row.names(pdx_degs)
        inter <- intersect(primary_degs, pdx_degs)
      }



      #  pdx
      degs_all = degs()$pdx

      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }


      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      degs_all$baseMeanNew <- 1 / (10^log((2^degs_all$AveExpr) + 1))

      if(input$keep_intersection == TRUE){ selectedGenes = inter }


      pv <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                             lab = rownames(degs_all),
                                             x = 'logFC',
                                             y = "baseMeanNew", #input$volc_y_axis,
                                             pCutoff =  1/10^as.numeric(input$basemean_cutoff),
                                             FCcutoff = timer()$lfc,
                                             xlab = bquote(~Log[2]~ 'fold change'),
                                             ylab = bquote(~Log[e]~ 'base mean + 1'),
                                             legendLabels = c('NS', expression(Log[2]~FC), 'Mean expr.', expression(Mean-expr~and~log[2]~FC)),
                                             title="MA plot -- PDX samples",
                                             selectLab = selectedGenes,
                                             pointSize=input$ma_pointSize,
                                             labSize=input$ma_labSize,
                                             colAlpha=input$ma_colAlpha,
                                             shape=as.numeric(input$ma_shape),
                                             legendPosition=input$ma_legendPosition,
                                             subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                             legendLabSize=input$ma_legendLabSize,
                                             legendIconSize=input$ma_legendIconSize,
                                             labFace=input$ma_labFace,
                                             boxedLabels=input$ma_boxedLabels,
                                             drawConnectors=input$ma_drawConnectors,
                                             gridlines.major = FALSE,
                                             gridlines.minor = FALSE,
                                             widthConnectors = 0.25,
                                             ylim = c(0, 10),
                                             col=input$ma_col) + coord_flip()


      #  primary
      degs_all = degs()$primary

      if(input$label_de_genes == TRUE){
        degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
        degs <- row.names(degs)
      } else {
        degs <- NULL
      }


      selectedGenes <- unique(c(input$multi_gene_symbol, degs))

      degs_all$baseMeanNew <- 1 / (10^log((2^degs_all$AveExpr) + 1))

      if(input$keep_intersection == TRUE){ selectedGenes = inter }

      pv2 <- EnhancedVolcano::EnhancedVolcano(degs_all,
                                              lab = rownames(degs_all),
                                              x = 'logFC',
                                              y = "baseMeanNew", #input$volc_y_axis,
                                              pCutoff =  1/10^as.numeric(input$basemean_cutoff),
                                              FCcutoff = timer()$lfc,
                                              xlab = bquote(~Log[2]~ 'fold change'),
                                              ylab = bquote(~Log[e]~ 'base mean + 1'),
                                              legendLabels = c('NS', expression(Log[2]~FC), 'Mean expr.', expression(Mean-expr~and~log[2]~FC)),
                                              title="MA plot -- primary samples",
                                              selectLab = selectedGenes,
                                              pointSize=input$ma_pointSize,
                                              labSize=input$ma_labSize,
                                              colAlpha=input$ma_colAlpha,
                                              shape=as.numeric(input$ma_shape),
                                              legendPosition=input$ma_legendPosition,
                                              subtitle=paste0(timer()$comp_a, " vs. ", timer()$comp_b),
                                              legendLabSize=input$ma_legendLabSize,
                                              legendIconSize=input$ma_legendIconSize,
                                              labFace=input$ma_labFace,
                                              boxedLabels=input$ma_boxedLabels,
                                              drawConnectors=input$ma_drawConnectors,
                                              gridlines.major = FALSE,
                                              gridlines.minor = FALSE,
                                              widthConnectors = 0.25,
                                              ylim = c(0, 10),
                                              col=input$ma_col) + coord_flip()
      grid.arrange(pv, pv2, nrow = 2)

    }


  }

  output$ma_plot <- renderCachedPlot({

      shiny::validate(
        need(make.names(timer()$comp_a) != make.names(timer()$comp_b),
             message = "Select two different groups to compare against.")
      )


      shiny::validate(
        need(length(input$ma_col) == 4,
             message = "Select exactly four colors corresponding to < abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff.")
			 )


      makeMaPlot()

  }, cacheKeyExpr = { list(degs(),
                           timer()$comp_a,
                           timer()$comp_b,
                           timer()$padj,
                           timer()$lfc,
                           timer()$top_n_from_de_genes,
                           input$label_de_genes,
                           input$multi_gene_symbol,
                           input$ma_drawConnectors,
                           input$ma_col,
                           input$ma_boxedLabels,
                           input$ma_labFace,
                           input$ma_legendIconSize,
                           input$ma_legendLabSize,
                           input$ma_legendPosition,
                           input$ma_shape,
                           input$ma_colAlpha,
                           input$ma_labSize,
                           input$ma_pointSize,
                           input$basemean_cutoff,
                           input$keep_intersection,
                           input$split_by_type)})


  output$deg_table <-  DT::renderDataTable({
    if(input$split_by_type == FALSE){
      degs_all = degs()
      degs <-  subset(degs_all, adj.P.Val < timer()$padj  & abs(logFC) > timer()$lfc) %>% as.data.frame()
      as.data.frame(degs)
    }
  },server=TRUE, rownames = TRUE, filter = list(position = 'top', clear = FALSE), options = list(scrollX = TRUE,
                   search=list(regex = TRUE, caseInsensitive = TRUE),
                   pageLength = 20))





  callModule(ggDownload, "ggvol", reactive({ makeVolPlot()}))
  callModule(ggDownload, "ggheat", reactive({makeHeatPlot()}))
  callModule(ggDownload, "ggma", reactive({makeMaPlot()}))

}


