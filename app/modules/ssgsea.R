
gseaInput <- function(id) {
  ## Define what's shown on the side for all plots of the Gene Exploration panel

  ns <- NS(id) # namespaced IDs; here: "ge" is a likely candidate

  title <- textInput(ns("title"),
                     label = "Plot title",
                     value = "")


  method <- selectizeInput(ns("method"),
                                "Select algorithm",
                                choices = c("gsva", "ssgsea"),
                                selected="gsva",
                                multiple=FALSE)


  gene_set_db <- selectizeInput(ns("gene_set_db"),
                                "Select MSigDB gene set collection (v6.2)",
                                choices =gsub(".gmt", "", list.files("msigdb/")),
                                selected=NULL,
                                multiple=FALSE)


  top_n <- numericInput(ns('top_n'), label = 'Include top n most variable gene sets: ', value = 25, min = 1,max = 2000, step=10)

 de_only <- checkboxInput(ns('de_only'),"Sort by statistical significance? Valid only when more than one subtype is selected.",value = TRUE)

   select_samples_by <- selectInput(inputId=ns("select_samples_by"),
                                   label="Select samples by",
                                   choices = NULL,
                                   selected = NULL)

  selected_samples <- checkboxGroupInput(inputId=ns("selected_samples"),
                                         label="",
                                         choices = NULL,
                                         selected = NULL)

#  remove_batch_effect <-selectInput(ns('remove_batch_effect'), label = "Remove batch effect?", choices = NULL, selected=NULL, multiple = FALSE)

select_all <- actionLink(ns("select_all"), "Select all", style = "color:#18cc7b")
deselect_all <- actionLink(ns("deselect_all"), "Deselect all", style = "color:#cc3c18")


  remove_samples <-selectInput(ns('remove_samples'), label = "Remove sample by name:", choices = NULL, selected=NULL, multiple = TRUE)

  include_samples <-selectInput(ns('include_samples'), label = "Include sample by name:", choices = NULL, selected=NULL, multiple = TRUE)


  # list of UIs to be exported to the add
  tagList(title,
          method,
          gene_set_db,
          top_n,
          de_only,
        #  hr(),
        #  remove_batch_effect,
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


gseaOutput <- function(id) {
  ns <- NS(id)
  tagList(
    shinydashboard::box( width = 12, title = "Gene set variation analysis heatmap",
                         fluidRow(column(width = 2,
                                         checkboxInput(ns('sample_labels_hm'),
                                                       "Display sample labels?",
                                                       value = TRUE)),
                                  column(width = 2,
                                         checkboxInput(ns('row_labels_hm'),
                                                       "Display pathway IDs?",
                                                       value = TRUE)),
                                  column(width = 2,
                                         checkboxInput(ns('trim_labels'),
                                                       "Trim pathway IDs?",
                                                       value = TRUE)),
                                  column(width=2,
                                         checkboxInput(ns('cluster_col'),
                                                       "Cluster columns?",
                                                       value = TRUE)),
                                  column(width=2,
                                         checkboxInput(ns('cluster_row'),
                                                       "Cluster rows?",
                                                       value = TRUE)),
                                  column(width = 2,
									  selectInput(ns('scale_by'),
									  label = 'Scale by:',
									  choices = c("none", "row", "column"),
									  selected = "none" ,multiple = FALSE))),
                         fluidRow(column(width = 3,
                                         selectInput(ns('column_annotation'), label = 'Annotate column by:', choices = NULL,
                                                     selected = NULL,multiple = TRUE)),
                                  column(width = 3,
                                         numericInput(ns('fontsize_row'), label = 'Row fontsize :', value = 8,
                                                      min = 1,step = 1)
                                  ),
                                  column(width=3,
                                  selectInput(ns('cutree_cols'), label = 'Cut cols into N groups:', choices = seq(1,20), selected = 1,
                                              multiple = FALSE)),
                                              column(width = 3,
                                                     selectInput(ns('color_palette_seq'),
                                                                 label = 'Select color palette:',
                                                                 choices = row.names(subset(RColorBrewer:::brewer.pal.info, category=="div"  | category=="seq")),
                                                                 selected = "RdYlBu" ,multiple = FALSE))),
                         plotOutput(ns("pheatmap_plot"), height="1200px"),
                         helpText("GSVA assesses the relative enrichment of gene sets across samples using a non-parametric approach. Conceptually, GSVA transforms a p-gene by n-sample gene expression matrix into a g-geneset by n-sample pathway enrichment matrix. This facilitates many forms of statistical analysis in the 'space' of pathways rather than genes, providing a higher level of interpretability.
                                  By default, the method is set to gsva (Hänzelmann et al, 2013) and the other option is ssgsea (Barbie et al, 2009)."),
									  ggDownloadUI(ns("ggheat")),
                         collapsible = FALSE))
}


gseaMod <- function(input, output, session, dds = dds, gsva = gsva, ssgsea=ssgsea, metadata) {

   updateSelectInput(session, "select_samples_by", choices = c(names(metadata)), selected="type")
   updateSelectInput(session, "facet_by_row", choices = c(none = ".", names(metadata)), selected=".")
   updateSelectInput(session, "facet_by_col", choices = c(none = ".", names(metadata)), selected=".")
   updateSelectInput(session, "include_samples", choices = c("",row.names(metadata)), selected=NULL)
   updateSelectInput(session, "column_annotation", choices = c("", names(metadata)), selected="type")
   updateSelectInput(session, "remove_batch_effect", choices = c("none","source", "batch", "patient", "purity"), selected=NULL)



   observeEvent(input$select_samples_by, {
     updateCheckboxGroupInput(session, "selected_samples", choices = unique(metadata[[input$select_samples_by]]), selected=input$selected_samples, inline=TRUE )
   }, ignoreInit = TRUE)




   mydds <- reactiveValues(retain_samples = NULL, remove_samples = NULL, include_samples = NULL, final_list_samples_to_include=NULL)
    timer <- reactive(list(final_list_samples_to_include = mydds$final_list_samples_to_include,
                           retain_samples  = mydds$retain_samples,
                           include_samples = mydds$include_samples)) %>% debounce(1000)

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



    heatmap_data <- reactive({
      shiny::validate(
        need(length(timer()$retain_samples) >0  | length(timer()$include_samples) > 0,
             message = "Select a sample to get started.")
      )

      withProgress(message = 'Extracting GSEA', {
        incProgress(0)

        selected_geneset <- input$gene_set_db

        if(input$method == "gsva"){
          df <- as.data.frame(gsva[[selected_geneset]])
        } else {
          df <- as.data.frame(ssgsea[[selected_geneset]])
        }
        df <- df[, timer()$final_list_samples_to_include]



        if(input$de_only == TRUE){
          if (length(levels(factor(colData(dds[,timer()$final_list_samples_to_include])$type))) > 1){
            design <- model.matrix(~factor(colData(dds[,timer()$final_list_samples_to_include])$type))
            fit <- lmFit(df, design)
            fit <- eBayes(fit)
            deGeneSets <- topTable(fit, coef=2:ncol(design), number=Inf, p.value = 1)
            df <- df[head(rownames(deGeneSets), input$top_n),]
          } else {
            df = df[head(order(rowVars(as.matrix(df)), decreasing = T), input$top_n),]
          }
        } else {
          df = df[head(order(rowVars(as.matrix(df)), decreasing = T), input$top_n),]
        }


        return(df)
      })
    }) %>%
      bindCache(timer()$final_list_samples_to_include,input$gene_set_db, input$method, input$de_only, input$top_n)




    makeHeatPlot <- function(){

      shiny::validate(
        need(length(timer()$final_list_samples_to_include) >0,
             message = "Select at least two samples to get started.")
      )

      heatmap_data = heatmap_data()

      if (length(input$column_annotation) > 0){
        expgroups <- as.data.frame(colData(dds)[,input$column_annotation])
        rownames(expgroups) <- colnames(dds)
        colnames(expgroups) <- input$column_annotation
      } else {
        expgroups = NULL
      }


      if(input$row_labels_hm == TRUE) {
        if(input$trim_labels == TRUE){
          labels_row=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_|^GO_|^HALLMARK|", "", rownames(heatmap_data))), 1, 40)
        } else {
          labels_row = row.names(heatmap_data)
        }
      }  else {labels_row = rep(" ", 10000) }


      pheatmap::pheatmap(heatmap_data,
                         cluster_rows = input$cluster_row,
                         cluster_cols = input$cluster_col,
                         cutree_cols = input$cutree_cols,
                         main=input$title,
                         show_colnames = input$sample_labels_hm,
                         labels_row=labels_row,
                         annotation_col = expgroups,
                         scale=input$scale_by,
                         color = colorRampPalette(rev(RColorBrewer:::brewer.pal(n = 7, name =input$color_palette_seq)))(100),
                         fontsize_row = input$fontsize_row)

    }


   output$pheatmap_plot <- renderCachedPlot({

       shiny::validate(
         need(length(timer()$final_list_samples_to_include) >0,
              message = "Select at least two samples to get started.")
       )


    makeHeatPlot()

         }, cacheKeyExpr = { list(heatmap_data(),
                 input$title,
                 input$cluster_row,
                 input$cluster_col,
				 input$cutree_cols,
				 input$sample_labels_hm,
				 input$scale_by,
				 input$fontsize_row,
				 input$trim_labels,
				 input$row_labels_hm,
				 input$column_annotation,
                 input$color_palette_seq,
                 timer()$final_list_samples_to_include,
                 input$gene_set_db,
                 input$method,
                 input$de_only,
                 input$top_n)})



   callModule(ggDownload, "ggheat", reactive({makeHeatPlot()}))
}









## gsub(".gmt", "", list.files("msigdb/"))
# gsva <- list()
# list.files("msigdb/")
# dds <- readRDS("./data/dds_march2019.RDS")
#vst_counts <- assay(dds, "vst")
#  gsva$biocarta <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/biocarta.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$canonical_pathways <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/canonical_pathways.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$chemical_and_genetic_perturbations <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/chemical_and_genetic_perturbations.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$go_biological_process <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_biological_process.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$go_cellular_component <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_cellular_component.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$go_molecular_function <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_molecular_function.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$hallmark <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/hallmark.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$immunologic_signatures <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/immunologic_signatures.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$kegg <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/kegg.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$oncogenic_signatures <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/oncogenic_signatures.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$reactome <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/reactome.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
#  gsva$transcription_factor_targets <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/transcription_factor_targets.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4)
# saveRDS(gsva, "gsva.RDS", compress=F)


# gsub(".gmt", "", list.files("msigdb/"))
# ssgsea <- list()
# list.files("msigdb/")
# dds <- readRDS("./data/dds_march2019.RDS")
# vst_counts <- assay(dds, "vst")
# ssgsea$biocarta <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/biocarta.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$canonical_pathways <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/canonical_pathways.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$chemical_and_genetic_perturbations <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/chemical_and_genetic_perturbations.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$go_biological_process <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_biological_process.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$go_cellular_component <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_cellular_component.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$go_molecular_function <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/go_molecular_function.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$hallmark <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/hallmark.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$immunologic_signatures <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/immunologic_signatures.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$kegg <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/kegg.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$oncogenic_signatures <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/oncogenic_signatures.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$reactome <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/reactome.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# ssgsea$transcription_factor_targets <- gsva(as.matrix(vst_counts), GSEABase::getGmt("msigdb/transcription_factor_targets.gmt"), min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, parallel.sz=4, method="ssgsea")
# saveRDS(ssgsea, "ssgsea.RDS", compress=F)
#
#
