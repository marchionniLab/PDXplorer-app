
diffInput <- function(id) {

  ns <- NS(id)


  title <- textInput(ns("title"), label = "Plot title", value = "")

  select_contrast_by <- selectInput(inputId=ns("select_contrast_by"),
                                   label="Select group",
                                   choices = NULL,
                                   selected = NULL)


  find_genes_by <-selectInput(ns('find_genes_by'), label = "Find DE genes by:", choices = c("contrast", "gene symbol"), selected="contrast", multiple = FALSE)

  selected_contrast <-selectInput(ns('selected_contrast'), label = "Select contrast:", choices = NULL, selected=NULL, multiple = FALSE)

  selected_multi_contrast <-selectInput(ns('selected_multi_contrast'), label = "Select contrasts:", choices = NULL, selected=NULL, multiple = TRUE)

  logfc_cutoff <-numericInput(ns('logfc_cutoff'), label = "Min. absolute log2-fold-change:", value = 0, min = 0, max = 1, step = 0.01)

  pval_cutoff <-numericInput(ns('pval_cutoff'), label = "Cutoff value for adj. p-values:", value = 0.05, min = 0.01, max = 20, step = 0.01)

  ontology <-selectInput(ns('ontology'), label = "Ontologies to include", choices = c("BP", "CC", "MF", "all"), selected="all", multiple = FALSE)


  gene_symbol <- selectizeInput(ns('gene_symbol'), label = "Find contrasts with following gene(s) DE", choices = NULL, selected=NULL, multiple = TRUE)



  tagList(select_contrast_by,
          find_genes_by,
          selected_contrast,
          gene_symbol,
          selected_multi_contrast,
          hr(),
          h5("Filter results"),
          logfc_cutoff,
          ontology,
          pval_cutoff

  )

}


diffOutput <- function(id) {
  ns <- NS(id)

  tagList(

    ## cluster
    tabsetPanel(id = ns("meta_tabs"),
                tabPanel("Differentially expressed genes",
    tabsetPanel(id = ns("tabs"),
                tabPanel("Scatter plot",
                         div(id = ns("scatter-outer"),
                         shinydashboard::box(title = "Scatter plot", width = 12, id = ns("scatterPlotBox"),
                                             plotOutput(ns("scatterPlot"), height="600px", click = clickOpts(id =ns("plot_click"))),
                                             collapsible = TRUE, collapsed = TRUE))),
                tabPanel("MD plot",
                         div(id = ns("md-outer"),
                         shinydashboard::box(title = "Mean-difference plot", width = 12, id = ns("MDplotBox"),
                                             plotOutput(ns("MDplot"), height="600px", click = clickOpts(id =ns("plot_click"))),
                                             collapsible = TRUE, collapsed = TRUE))),
                tabPanel("Volcano plot",
                         div(id = ns("volcano-outer"),
                         shinydashboard::box(title = "Volcano plot", width = 12, id=ns("volcanoPlotBox"),
                                             plotOutput(ns("volcanoPlot"), height="600px", click = clickOpts(id =ns("plot_click"))),
                                             collapsible = TRUE, collapsed = TRUE)))
    ),
  #  actionLink(ns("copy_to_clipboard"), "Copy gene list to clipboard"),
    uiOutput(ns("weblink_to_heatmap")),
    hr(),
    DT::dataTableOutput(ns('deg_table'), height="1000px")

# helpText("If a gene is detected as differentially expressed in any of the contrasts with the adjusted p.value and log-fold change specified to the left, each contrast name will appear here along with corresponding log-fold-change values.")
                ),
  tabPanel("Intersections",
           plotOutput(ns("upsetr"), height="600px"),
           DT::dataTableOutput(ns('upsetr_table'), height="1000px")
  ),
  tabPanel("GO enrichments",
           DT::dataTableOutput(ns('go_table'),  height="1000px")
  ),
  tabPanel("Kegg pathway analysis",
           shinydashboard::box(title = "Barcode plot", width = 12,
                               plotOutput(ns("barcodePlot"),height="350px", click = clickOpts(id =ns("plot_click"))),
                               verbatimTextOutput(ns("barcode_caption")),
                               collapsible = TRUE, collapsed = TRUE),
           DT::dataTableOutput(ns('kegg_table'),  height="1000px")
           #         DT::dataTableOutput(ns('go_table'),  height="1000px")
  )
)
    )
}


diffMod <- function(input, output, session, dds, metadata, deg, go, keggNames, keggLinks, ann) {

#  row.names(deg$primary$p.value[apply(deg$primary$p.value[, , drop=F]<0.05, 1, any),])
  # row.names(deg[[input$select_contrast_by]]$p.value)



  updateSelectInput(session, "select_contrast_by", choices=names(deg)[!(names(deg) %in% "decoder")], selected="primary")

  updateSelectInput(session, "selected_contrast", choices=deg$decoder$real_compare[deg$decoder$valid_compare %in% c(colnames(deg[["primary"]]$coefficients))])

  updateSelectInput(session, "selected_multi_contrast", choices=c(paste0("PDX_" , colnames(deg[["PDX"]]$coefficients)), paste0("primary_" , colnames(deg[["primary"]]$coefficients))))

 # updateSelectizeInput(session, "gene_symbol", choices=row.names(deg[["primary"]]$p.value[apply(deg[["primary"]]$p.value[, , drop=F]<0.05, 1, any),]), server=TRUE)

  shinyjs::hide("gene_symbol")

  observeEvent(input$select_contrast_by, {
    updateSelectInput(session, "selected_contrast", choices=deg$decoder$real_compare[deg$decoder$valid_compare %in% c(colnames(deg[[input$select_contrast_by]]$coefficients))], selected=input$selected_contrast)
      updateSelectizeInput(session, "gene_symbol", choices=row.names(deg[[input$select_contrast_by]]$p.value[apply(deg[[input$select_contrast_by]]$p.value[, , drop=F]<0.05, 1, any),]), server=TRUE)

 }, ignoreInit = TRUE)

  observeEvent(input$find_genes_by, {
    if (input$find_genes_by == "gene symbol"){
      shinyjs::show("gene_symbol")
      shinyjs::hide("selected_contrast")
      shinyjs::hide("tabs")
      shinyjs::hide("weblink_to_heatmap")
      shinyjs::hide("scatterPlot")
      shinyjs::hide("MDplot")
      shinyjs::hide("volcanoPlot")
      shinyjs::hide("scatterPlotBox")
      shinyjs::hide("MDplotBox")
      shinyjs::hide("volcanoPlotBox")
      shinyjs::hide("volcano-outer")
      shinyjs::hide("md-outer")
      shinyjs::hide("scatter-outer")
    #  shinyjs::show("helptxtdiv")
      } else {
      shinyjs::hide("gene_symbol")
      shinyjs::show("selected_contrast")
      shinyjs::show("tabs")
      shinyjs::show("weblink_to_heatmap")
      shinyjs::show("scatterPlot")
      shinyjs::show("MDplot")
      shinyjs::show("volcanoPlot")
  #    shinyjs::show("scatterPlotBox")
  #    shinyjs::show("MDplotBox")
  #    shinyjs::show("volcanoPlotBox")
      shinyjs::show("volcano-outer")
      shinyjs::show("md-outer")
      shinyjs::show("scatter-outer")
     # shinyjs::hide("helptxtdiv")
    }
  }, ignoreInit = TRUE)


  data <- reactiveValues(deg = NULL, kegg=NULL, inter = NULL)

  observeEvent(input$meta_tabs, {
    if (input$meta_tabs ==  "Intersections"){
      shinyjs::hide("selected_contrast")
      shinyjs::hide("select_contrast_by")
      shinyjs::show("selected_multi_contrast")
      shinyjs::hide("ontology")
      shinyjs::show("logfc_cutoff")
      shinyjs::hide("gene_symbol")
      shinyjs::hide("find_genes_by")
      } else if (input$meta_tabs ==  "Differentially expressed genes") {
        shinyjs::hide("selected_multi_contrast")
       # shinyjs::show("selected_contrast")
        shinyjs::show("select_contrast_by")
        shinyjs::show("find_genes_by")

        shinyjs::hide("ontology")
        shinyjs::show("logfc_cutoff")
     #   shinyjs::show("gene_symbol")

        if (input$find_genes_by == "gene symbol"){
          shinyjs::show("gene_symbol")
          shinyjs::hide("selected_contrast")
          shinyjs::hide("tabs")
          shinyjs::hide("weblink_to_heatmap")
          shinyjs::hide("scatterPlot")
          shinyjs::hide("MDplot")
          shinyjs::hide("volcanoPlot")
          shinyjs::hide("scatterPlotBox")
          shinyjs::hide("MDplotBox")
          shinyjs::hide("volcanoPlotBox")
          } else {
          shinyjs::hide("gene_symbol")
          shinyjs::show("selected_contrast")
          shinyjs::show("tabs")
          shinyjs::show("weblink_to_heatmap")
          shinyjs::show("scatterPlot")
          shinyjs::show("MDplot")
          shinyjs::show("volcanoPlot")
      #    shinyjs::show("scatterPlotBox")
      #    shinyjs::show("MDplotBox")
       #   shinyjs::show("volcanoPlotBox")
        }
      } else if (input$meta_tabs ==  "GO enrichments") {
        shinyjs::show("ontology")
        shinyjs::hide("selected_multi_contrast")
        shinyjs::show("selected_contrast")
        shinyjs::show("select_contrast_by")
        shinyjs::hide("logfc_cutoff")
        shinyjs::hide("gene_symbol")
        shinyjs::hide("find_genes_by")

      } else if (input$meta_tabs == "Kegg pathway analysis"){
        shinyjs::hide("logfc_cutoff")
        shinyjs::hide("ontology")
        shinyjs::hide("selected_multi_contrast")
        shinyjs::show("selected_contrast")
        shinyjs::show("select_contrast_by")
        shinyjs::hide("gene_symbol")
        shinyjs::hide("find_genes_by")

      }
  })

  output$upsetr <- renderPlot({
    shiny::validate(
      need(length(input$selected_multi_contrast) > 1 ,
           message = "Select at leasts two contrasts to get started.")
    )
    # https://github.com/hms-dbmi/UpSetR/issues/85

    fromList <- function (input) {
      # Same as original fromList()...
      elements <- unique(unlist(input))
      data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
      }))
      data[is.na(data)] <- as.integer(0)
      data[data != 0] <- as.integer(1)
      data <- data.frame(matrix(data, ncol = length(input), byrow = F))
      data <- data[which(rowSums(data) != 0), ]
      names(data) <- names(input)
      # ... Except now it conserves your original value names!
      row.names(data) <- elements
      return(data)
    }
    res <- list()
    res = sapply(input$selected_multi_contrast, function(x) row.names( topTable(deg[[strsplit(x, "_")[[1]][1] ]], coef=strsplit(x, "_")[[1]][2], adjust="BH", p.value=input$pval_cutoff, lfc=input$logfc_cutoff, n=Inf)))
    UpSetR::upset(fromList(res), order.by = "freq", nsets=50,  main.bar.color="#E2463E",sets.bar.color="#26699A", text.scale=1)

  data$inter <- row.names(fromList(res))

  })

  output$upsetr_table <- DT::renderDataTable({

    shiny::validate(
      need(length(input$selected_multi_contrast) > 1 ,
           message = "")
    )

    res = lapply(input$selected_multi_contrast, function(x)  topTable(deg[[strsplit(x, "_")[[1]][1] ]], coef=strsplit(x, "_")[[1]][2], adjust="BH", p.value=input$pval_cutoff, lfc=input$logfc_cutoff, n=Inf)[,c(1,3,5)])
    names(res) <- input$selected_multi_contrast
    df <- res %>%  Reduce(function(dtf1,dtf2) inner_join(dtf1,dtf2,by="SYMBOL"), .)
    cleaned_df <- df[,c(1,2, grep("logFC", names(df)))]
    colnames(cleaned_df)[1:2] <- c("Gene", "Description")
    colnames(cleaned_df)[-c(1,2)] <- paste0("logFC.", input$selected_multi_contrast)
    cleaned_df

    }, escape = FALSE, selection = list(mode = 'single', target = 'row'), server=TRUE, rownames = FALSE, options = list(scrollX = TRUE, search=list(regex = TRUE, caseInsensitive = TRUE),pageLength = 20,selection = list(mode = 'single', target = 'cell')) )



  observe({
#    clipr::write_clip(row.names(data$deg))
    shiny::validate(
      need(!is.null(data$deg),
           message = "Select a contrast to get started.")
    )
    df <- topTable(deg[[input$select_contrast_by]], coef=as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare), adjust="BH",  p.value=input$pval_cutoff, lfc=input$logfc_cutoff, n=50)
    df$Gene <- row.names(df)
    df <- df[,c(7,1,2,5)]

    first_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][1]
    second_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][2]
    first_sel <- gsub("\\+$",  "%2B", first_sel)
    second_sel <- gsub("\\+$",  "%2B", second_sel)
    type=paste0("&ge-select_samples_by=%22tumor%22&ge-selected_samples=",paste0("%5B%22",first_sel,"%22%2C%22",second_sel,"%22%5D&"))

    first_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][1]
    second_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][2]
    sel_sample = paste0("&ge-remove_samples=%5B%22", paste( row.names(  subset(metadata,tumor %in% c(first_sel, second_sel) & type != input$select_contrast_by) ),"%22%2C%22",collapse="", sep="") , "%22%5D&")

    gene=paste0("&ge-multi_gene_symbol=%5B%22", paste0(row.names(df), collapse="%22%2C%22"),"%22%5D")

#gene=paste0("&ge-gene_symbol=%22",df$Gene,"%22")
    urls <- paste0("https://abc.med.cornell.edu/shiny/PDXplorer-rna/?_inputs_&ge-color_palette=%22npg%22&ge-facet_by_col=%22.%22&ge-facet_by_row=%22.%22&ge-group_by=%5B%22type%22%2C%22tumor%22%5D&ge-include_samples=null&ge-legend_position=%22right%22&ge-plot_click_ge=null&ge-plot_style=%22boxplot%22&ge-remove_samples=null&ge-sample_labels=false&ge-select_all=0&ge-select_samples_by=%22type%22&ge-title=%22%22&ge-ylimZero=true&tabs=%22Gene%20exploration%22", type, sel_sample, gene, "&ge-tabs=%22Heatmap%22&ge-column_annotation=%22tumor%22")
    refs <- paste0(paste0(paste0("<a href='",  urls, "' target='_blank'>"), "", "Send top 50 genes to heatmap"), "", "</a>")
    output$weblink_to_heatmap <- renderUI({
      HTML(refs)
    })
  })


  output$deg_table <-  DT::renderDataTable({

    shiny::validate(
      need(input$selected_contrast != "",
           message = "")
    )

if (input$find_genes_by == "contrast"){
  ## full df
  df <- topTable(deg[[input$select_contrast_by]], coef=as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare), adjust="BH",  p.value=1, lfc=0, n=Inf)
  df$Gene <- row.names(df)
  # df <- df[,c(7,1,2,5)]
  data$deg <- df


  df <- topTable(deg[[input$select_contrast_by]], coef=as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare), adjust="BH",  p.value=input$pval_cutoff, lfc=input$logfc_cutoff, n=Inf)
  df$Gene <- row.names(df)
  df <- df[,c(11,3,5,6,9)]
  colnames(df)[2] <- "Description"

  first_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][1]
  second_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][2]
  first_sel <- gsub("\\+$",  "%2B", first_sel)
  second_sel <- gsub("\\+$",  "%2B", second_sel)
  type=paste0("&ge-select_samples_by=%22tumor%22&ge-selected_samples=",paste0("%5B%22",first_sel,"%22%2C%22",second_sel,"%22%5D&"))

  first_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][1]
  second_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][2]
  sel_sample = paste0("&ge-remove_samples=%5B%22", paste( row.names(  subset(metadata,tumor %in% c(first_sel, second_sel) & type != input$select_contrast_by) ),"%22%2C%22",collapse="", sep="") , "%22%5D&")
  gene=paste0("&ge-gene_symbol=%22",df$Gene,"%22")
  urls <- paste0("https://abc.med.cornell.edu/shiny/PDXplorer-rna/?_inputs_&ge-color_palette=%22npg%22&ge-facet_by_col=%22.%22&ge-facet_by_row=%22.%22&ge-group_by=%5B%22type%22%2C%22tumor%22%5D&ge-include_samples=null&ge-legend_position=%22right%22&ge-plot_click_ge=null&ge-plot_style=%22boxplot%22&ge-remove_samples=null&ge-sample_labels=false&ge-select_all=0&ge-select_samples_by=%22type%22&ge-title=%22%22&ge-ylimZero=true&tabs=%22Gene%20exploration%22", type, sel_sample, gene)
  refs <- paste0(paste0(paste0("<a href='",  urls, "' target='_blank'>"), "", df$Gene), "", "</a>")
  df$Gene <- refs
  as.data.frame(df)
} else {


  shiny::validate(
    need(input$gene_symbol != "",
         message = "Select a gene(s) to get started.")
  )

  df <- lapply(input$gene_symbol, function(x) {
    col_names <- colnames(deg[[input$select_contrast_by]]$p.value[x, , drop=F])[deg[[input$select_contrast_by]]$p.value[x, , drop=F] < input$pval_cutoff]
    logfc_names <-  colnames(deg[[input$select_contrast_by]]$coefficients[x,col_names, drop=F])[abs(deg[[input$select_contrast_by]]$coefficients[x,col_names, drop=F]) > input$logfc_cutoff]
    as.data.frame(t(deg[[input$select_contrast_by]]$coefficients[x,logfc_names]))
  }
  ) %>% plyr::rbind.fill()
  row.names(df) <- input$gene_symbol
  df <- t(df)

}

  }, escape = FALSE, selection = list(mode = 'single', target = 'row'), server=TRUE, rownames = TRUE, options = list(scrollX = TRUE, search=list(regex = TRUE, caseInsensitive = TRUE),pageLength = 20,selection = list(mode = 'single', target = 'cell')) )


  output$volcanoPlot <- renderPlot({
    shiny::validate(
      need(!is.null(data$deg),
           message = "Select a contrast to get started.")
    )

    gene_list <- data$deg
    gene_list$threshold = as.factor(abs(gene_list$logFC) > input$logfc_cutoff & gene_list$adj.P.Val < input$pval_cutoff)
    p <- ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold, key=row.names(gene_list))) +
      geom_point(size=0.4) +
      xlim(c(min(gene_list$logFC), max(gene_list$logFC))) +
      ylab("-log10(FDR)") +
      xlab("log2(fold change)")  + theme_bw(base_size=14) +
      theme(legend.position="bottom") +
      scale_colour_manual(name = "match filters", values = c("#808080", "#D75F1C")) + geom_hline(yintercept = -log10(as.numeric(input$pval_cutoff)), linetype="dotted") + geom_vline(aes(xintercept=-as.numeric((input$logfc_cutoff))), linetype="solid") + geom_vline(aes(xintercept=input$logfc_cutoff), linetype="dashed") + geom_vline(aes(xintercept=0), linetype="solid", size=0.2, show.legend=FALSE)  + ggtitle(input$selected_contrast)
    p
  })


  output$MDplot <- renderPlot({
    shiny::validate(
      need(!is.null(data$deg),
           message = "Select a contrast to get started.")
    )

    gene_list <- data$deg
    gene_list$threshold = as.factor(abs(gene_list$logFC) > input$logfc_cutoff & gene_list$adj.P.Val < input$pval_cutoff)
    p <- ggplot(data=gene_list, aes(x=AveExpr, y=logFC, colour=threshold, key=row.names(gene_list))) +
      geom_point(size=0.4) +
      ylim(c(min(gene_list$logFC), max(gene_list$logFC))) +
      xlab("Mean expression") +
      ylab("log2(fold change)") +   theme_bw(base_size=14) +
      theme(legend.position="bottom") +
      scale_colour_manual(name = "match filters", values = c("#808080", "#D75F1C"))  + geom_hline(aes(yintercept=0), linetype="solid", size=0.2) + geom_vline(aes(xintercept=0), linetype="solid", size=0.2) + geom_hline(aes(yintercept=-as.numeric((input$logfc_cutoff))), linetype="solid") + geom_hline(aes(yintercept=input$logfc_cutoff), linetype="dashed")  + ggtitle(input$selected_contrast)
    p
  })



  output$scatterPlot <- renderPlot({
    shiny::validate(
      need(!is.null(data$deg),
           message = "Select a contrast to get started.")
    )

    first_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][1]
    second_sel <- strsplit(as.character(subset(deg$decoder, real_compare == input$selected_contrast)$real_compare), " vs. ")[[1]][2]

    if(length(row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(first_sel))) ) > 1){
      first_mean <- log2(rowMeans(assay(dds)[,  row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(first_sel))) ]) + 0.5)
    } else {
      first_mean <- log2(assay(dds)[, row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(first_sel))) ] + 0.5)
    }

    if(length(row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(second_sel)))) > 1){
      second_mean <- log2(rowMeans(assay(dds)[,  row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(second_sel))) ]) + 0.5)
    } else {
      second_mean <- log2(assay(dds)[, row.names(subset(metadata, type == input$select_contrast_by & tumor %in% c(second_sel))) ] + 0.5)
    }

    df <- data.frame(first_sel = first_mean, second_sel = second_mean)
    names(df) <- c(first_sel, second_sel)
    row.names(df) <- row.names(assay(dds))
    df <- df[row.names(deg$primary),]

    gene_list <- data$deg
    gene_list <- gene_list[row.names(df), ]
    threshold = as.factor(abs(gene_list$logFC) > input$logfc_cutoff & gene_list$adj.P.Val < input$pval_cutoff)

    p <- ggplot(data=df,  aes(x=get(first_sel), y=get(second_sel), colour=threshold)) +
      geom_point(size=0.4) +
    #  ylim(c(min(gene_list$logFC), max(gene_list$logFC))) +
      xlab( paste0("log2(AveExpr) - ", first_sel)) +
      ylab( paste0("log2(AveExpr) - ", second_sel)) +   theme_bw(base_size=14) +
      theme(legend.position="bottom") +
      scale_colour_manual(name = "match filters", values = c("#808080", "#D75F1C"))  + geom_hline(aes(yintercept=0), linetype="solid", size=0.2) + geom_vline(aes(xintercept=0), linetype="solid", size=0.2)  + geom_abline(intercept=as.numeric(input$logfc_cutoff), linetype="dashed") + geom_abline(intercept=-as.numeric(input$logfc_cutoff), linetype="solid") + geom_abline(intercept=0, linetype="solid") + ggtitle(input$selected_contrast)
    p
  })


  output$go_table <-  DT::renderDataTable({

    shiny::validate(
      need(input$selected_contrast != "",
           message = "")
    )
    if(input$ontology == "all"){ selected_ontology = c("BP", "CC", "MF") } else {selected_ontology = input$ontology}

    subset(topGO(go[["go"]][[input$select_contrast_by]][[as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare)]], ontology=selected_ontology, n=Inf), P.Up < input$pval_cutoff | P.Down < input$pval_cutoff)

  }, escape = FALSE, selection = list(mode = 'single', target = 'row'), server=TRUE, rownames = FALSE, options = list(scrollX = TRUE, search=list(regex = TRUE, caseInsensitive = TRUE),pageLength = 20,selection = list(mode = 'single', target = 'cell')) )




  output$kegg_table <-  DT::renderDataTable({

    shiny::validate(
      need(input$selected_contrast != "",
           message = "")
    )
    if(input$ontology == "all"){ selected_ontology = c("BP", "CC", "MF") } else {selected_ontology = input$ontology}

    df <- subset(topKEGG(go[["kegg"]][[input$select_contrast_by]][[as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare)]], n=Inf), P.Up < input$pval_cutoff | P.Down < input$pval_cutoff)
    data$kegg <- df

    df

  }, escape = FALSE, selection = list(mode = 'single', target = 'row'), server=TRUE, rownames = FALSE, options = list(scrollX = TRUE, search=list(regex = TRUE, caseInsensitive = TRUE),pageLength = 20,selection = list(mode = 'single', target = 'cell')) )


  output$barcodePlot <- renderPlot({
    shiny::validate(
      need(!is.null(data$kegg),
           message = "Select a contrast to get started.")
    )

    index <- row.names(deg[[input$select_contrast_by]]) %in% ann[ann$ENTREZID %in% subset(keggLinks, PathwayID == row.names(data$kegg[input$kegg_table_rows_selected,]))$GeneID, ]$SYMBOL
    barcodeplot(deg[[input$select_contrast_by]]$coeff[,as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare)], index=index, main=paste0("LogFC: ", subset(keggNames, PathwayID == row.names(data$kegg[input$kegg_table_rows_selected,]))$Description), labels=c(strsplit( input$selected_contrast, " vs. ")[[1]][1], strsplit( input$selected_contrast, " vs. ")[[1]][2]))

  })

  output$barcode_caption <- renderPrint({ "Select a row in the table below to get started.  In the plot, all genes are ranked from left to right by decreasing log-fold change for the contrast and the genes within the gene set are represented by vertical bars, forming the barcode-like pattern. The curve (or worm) above the barcode shows the relative local enrichment of the bars in each part of the plot. The dotted horizontal line indicates neutral enrichment; the worm above the dotted line shows enrichment while the worm below the dotted line shows depletion." })

  # observeEvent(input$kegg_table_rows_selected, {
  #   print(input$kegg_table_rows_selected)
  #   print( row.names(data$kegg[input$kegg_table_rows_selected,]) )
  #   index <- row.names(deg[[input$select_contrast_by]]) %in% ann[ann$ENTREZID %in% subset(keggLinks, PathwayID == row.names(data$kegg[input$kegg_table_rows_selected,]))$GeneID, ]$SYMBOL
  #   print(table(index))
  #   barcodeplot(deg[[input$select_contrast_by]]$coeff[,as.character(subset(deg$decoder, real_compare == input$selected_contrast)$valid_compare)], index=index, main=paste0("LogFC: ", subset(keggNames, PathwayID == row.names(data$kegg[input$kegg_table_rows_selected,]))$Description), labels=c(" "," "))
  #
  # })

}




