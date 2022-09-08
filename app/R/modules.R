## additional argument to define the label on the downloadButton
pdfDownloadUI <- function(id, label = "Download PDF") {
  ns <- shiny::NS(id)
  downloadButton(ns("download"), label)
}

## allow users of the module to input a (reactive) data.frame to download as csv and a name for the file
pdfDownload <- function(input, output, session, makePlot,  plotHeight = plotHeight, plotWidth = plotWidth,
                        filename = paste0("plot_", format(Sys.time(), "%a-%b-%d-%H:%M:%S-%Y"), ".pdf")) {
  
  if (is.reactive(plotHeight)) {
    plotHeight <- as.numeric(plotHeight())
  }
  
  if (is.reactive(plotWidth)) {
    plotWidth = as.numeric(plotWidth())
  }
  
  
  output$download <- downloadHandler(
    filename = function() {
      filename
    },
    content = function(file) {
   #   pdf(file)  #    dev.off()
      
   ggsave(filename = file, plot=makePlot(), device="pdf", width=plotWidth, height=plotHeight)
    }
  )
}