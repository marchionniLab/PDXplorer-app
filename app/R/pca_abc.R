## for saving the coordinates of PCA or tSNE, see export_import.R
# ...............................................................................


#' Prepare the expression matrix for PCA analysis
#'
#' @param object a matrix of values where rows are genes (features) and columns
#' are samples.
#' @param n_VarGenes numeric scalar indicating the number of the most
#' variable features to use for the PCA.
#' @param scale_features logical, should the expression values be standardised
#' so that each feature has unit variance? Default is \code{TRUE}.
#'
#' @return a matrix where rows (x) are cells and columns (y) are genes
#' @examples
#' \dontrun{
#' pca_mat <- prep_PCA(assay(myrlt),
#'   n_VarGenes = input$pca_nrgenes,
#'   scale_features = TRUE
#' )
#' pcres <- prcomp(pca_mat)
#' pccoords <- return_PCA_coords(pcres,
#'   cell_names = colnames(myrlt),
#'   n_PCs = 20
#' )
#' }
prep_PCA <- function(object, n_VarGenes, scale_features) {
  exprs_mat <- object

  ## Define features to use: either n_VarGenes, or if a set of features is defined,
  ## then those
  rv <- matrixStats::rowVars(exprs_mat)
  feature_set <- order(rv, decreasing = TRUE)[seq_len(min(n_VarGenes, length(rv)))]

  # Subsetting to the desired features (do NOT move below 'scale()')
  exprs_for_pca <- exprs_mat[feature_set, , drop = FALSE]

  ## Standardise expression if scale_features argument is TRUE
  if (scale_features) {
    exprs_for_pca <- scale(t(exprs_for_pca))
  }

  ## Drop any features with zero variance
  keep_feature <- (matrixStats::colVars(exprs_for_pca) > 0.001)
  keep_feature[is.na(keep_feature)] <- FALSE
  exprs_for_pca <- exprs_for_pca[, keep_feature]

  return(exprs_for_pca)
}


#' Return PCA coordinates
#'
#' @description This function takes the output of \code{prcomp()} and
#' returns a list with two data.frames, one for the loadings of the genes,
#' one for the values of the eigenvector multiplication (cells).
#'
#' @details This function assumes that information about cells were
#' originally stored in rows and information about genes were kept
#' in columns.
#'
#' The \code{attr} slot of the data.frame for the cells will contain \code{
#' 	"percentVar"} and \code{"variance"}, which can be used for additional
#' diagnostic plots.
#'
#' @param pca_results result of \code{prcomp()}
#' @param cell_names vector of names for the rows of the matrix supplied
#' to \code{prcomp()}. Can be NULL (default).
#' @param n_PCs number of PCs to be returned; can be set to \code{Inf} to
#' get all (not recommended since this will equal the number of rows of the
#' 	original matrix)
#' @return list with two data.frames, one for the loadings of the genes,
#' one for the values of the eigenvector multiplication (cells)
#' @examples
#' \dontrun{
#' pca_mat <- prep_PCA(assay(myrlt),
#'   n_VarGenes = input$pca_nrgenes,
#'   scale_features = TRUE
#' )
#' pcres <- prcomp(pca_mat)
#' pccoords <- return_PCA_coords(pcres,
#'   cell_names = colnames(myrlt),
#'   n_PCs = 20
#' )
#' }
return_PCA_coords <- function(pca_results, cell_names = NULL, n_PCs) {
  if (class(pca_results) != "prcomp") {
    stop("The pca results are not of type 'prcomp'. Rerun prcomp().")
  }

  # define the number of PCs that are going to be returned
  ncomp_out <- min(n_PCs, ncol(pca_results$x))

  # extract info about variance captured by the PCs
  percentVar <- pca_results$sdev^2 / sum(pca_results$sdev^2)
  variance <- pca_results$sdev^2

  # cells
  df_out <- pca_results$x[, 1:ncomp_out]
  if (!is.null(cell_names)) {
    rownames(df_out) <- cell_names
  }
  attr(df_out, "percentVar") <- percentVar[1:ncomp_out]
  attr(df_out, "variance") <- variance[1:ncomp_out]

  # genes
  df_out2 <- data.frame(pca_results$rotation[, 1:ncomp_out])

  return(list(cells = df_out, genes = df_out2))
}



## Dealing with PCA results-----------------------------------------------------

#' Dot plot for the loadings of the top contributors to individual PCs
#'
#' @description Similar in spirit to Seurat::VizPlot(); wrapper function for
#' \code{\link{get_PCA_contributors}} and \code{\link{make_PCA_contributor_plot}}.
#'
#' @param in_df data.frame of PC loadings for genes or cells
#' @param n_pc number of PCs to plot
#' @param n_factors number of genes/cells to plot; the top \code{n_factors}
#' genes/cells for eac PC will be shown
#' @param return_df logical indicating whether the data.frame used for plotting
#' should be returned (default: FALSE)
#'
#' @return ggplot (and data.table if wanted)
#' @seealso \code{\link{Seurat::VizPCA}}, \code{\link{make_PCA_ABC}}
#' @examples
#' \dontrun{
#' plot_PCA_contributors(pca_medScaled$genes, n_pc = 5, n_factors = 10)
#' }
#'
#' @export
#'
plot_PCA_contributors <- function(in_df, n_pc = 1, n_factors = 10, return_df = FALSE) {
  top_contributors <- get_PCA_contributors(in_df, n_pc, n_factors)

  P <- make_PCA_contributor_plot(top_contributors)
  return(P)

  if (return_df) {
    return(top_contributors)
  }
}

#' Determine the top contributors for individual PC loadings
#'
#' @description Wrapper function for \code{\link{get_top_contributor}} to
#' get the top contributors to multiple PC loadings in a skinny data.frame
#' format that's fit for ggplot2.
#'
#' @param in_df data.frame of PC loadings for genes or cells
#' @param n_pc number of PC(s) for which to extract the top contributor(s)
#' @param n_factors how many top contributors should be extracted per PC
#'
#' @return data.frame in long format with "id", "value", "PC"
#' @examples
#' \dontrun{
#' top_contributors <- get_PCA_contributors(in_df, n_pc, n_factors)
#' P <- make_PCA_contributor_plot(top_contributors)
#' print(P)
#' }
#' #' @seealso \code{\link{get_top_contributor}}, \code{\link{make_PCA_contributor_plot}}
get_PCA_contributors <- function(in_df, n_pc = 1, n_factors = 10) {
  # define the PCs
  pcs <- names(in_df)[c(1:n_pc)]

  # get a long data.frame
  pca_df <- get_long_PCdf(in_df)

  # get the top contributors for each PC
  tops <- lapply(pcs, function(x) {
    get_top_contributor(pca_df,
      which_PC = x,
      n_out = n_factors,
      return_df = TRUE
    )
  })
  # tops <- data.table::rbindlist(tops)
  tops <- do.call("rbind", tops)

  return(tops)
}

#' Generate point plots for PCA contributors
#'
#' @description Short wrapper around the ggplot2 routine for making a dot plot
#' of PCA contributors and their loadings.
#'
#' @param tops long data.frame with the following columns: "id", "value", "PC"
#' @return ggplot2 object
#' @seealso \code{\link{plot_PCA_contributors}}, \code{\link{get_PCA_contributors}}
#' @examples
#' \dontrun{
#' tops <- get_PCA_contributors(in_df, n_pc, n_factors)
#' contri_plot <- make_PCA_contributor_plot(tops)
#' print(contri_plot)
#' }
make_PCA_contributor_plot <- function(tops) {
  ABCutilities::check_columns(c("id", "value", "PC"), tops, "top contributor df", "make_PCA_contributor_plot")

  min_load <- min(tops$value)
  max_load <- max(tops$value)

  tops$id <- reorder(tops$id, tops$value)

  P <- ggplot(tops, aes(x = id, y = value)) +
    geom_point(size = 3) +
    coord_cartesian(ylim = c(min_load, max_load)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~PC, scales = "free") +
    ylab("PC loading") +
    xlab("")

  return(P)
}


#' Helper function: turning the PCA data.frame into a long one
#'
#' @param pca_df data.frame with the component loadings for the PCs for
#' either genes or cells; rownames should have the gene or cell names, colnames
#' correspond to the individual PCs
#' @return data.frame in skinny format (one row per gene|cell ~ loading combi)
#' @examples
#' \dontrun{
#' get_long_PCdf(pca_result$genes)
#' }
get_long_PCdf <- function(pca_df) {
  # pc <- as.data.table(pca_df)
  pc <- pca_df
  attributes(pc)[names(attributes(pc)) %in% c("percentVar", "variance")] <- NULL
  pc <- data.frame(pc)
  pc$id <- rownames(pca_df)
  # pc <- melt.data.table(pc, id.vars = "id", variable.name = "PC")
  pc <- reshape2::melt(pc, id.vars = "id", variable.name = "PC")
  pc$PC <- as.character(pc$PC)
  return(pc)
}


#' Extract the genes or cells that contribute most to an individual principal
#' component
#'
#' @description component loading
#'
#' @param pca_df data.frame in skinny format with columns "id", "PC", "value"
#' @param which_PC indicate the name of the PC for which the top contributor(s)
#' should be returned (default: "PC1")
#' @param n_out number of genes/cell to be returned (default: 6)
#' @param return_dt logical indicating whether the values should be returned in
#' addition to the names (default: FALSE)
#'
#' @return vector of cell/gene names
get_top_contributor <- function(pca_df, which_PC = "PC1", n_out = 6, return_df = FALSE) {
  if (!"PC" %in% names(pca_df)) {
    stop(
      "The df you supplied does not contain a column named 'PC'. ",
      "Make sure get_PCA_contributors() instead which takes care of that."
    )
  }
  p_out <- subset(pca_df, PC == which_PC)
  p_out$abs_val <- abs(p_out$value)
  # print(head(p_out))
  # print(str(p_out))
  p_out <- p_out[order(-p_out$abs_val), ]

  if (!return_df) {
    out <- head(p_out$id, n = n_out)
  } else {
    setkey(p_out, NULL)
    out <- head(p_out, n = n_out)
  }
  return(out)
}

#' Principal components (cor)relation with experimental covariates
#'
#' @description Based on \code{pcaexplorer::correlatePCs}. Computes the
#' significance of (cor)relations between PCA scores and the sample
#' experimental covariates.
#'
#' @details Based on \code{pcaexplorer::correlatePCs}.
#' For comparing the PCA scores with \emph{categorical variables},
#' the Kruskal-Wallis test ist used.
#' For \emph{continuous variables}, the \code{cor.test} based on Spearman's
#' correlation is used.
#'
#' @param pca_res_samples The PCA result for the cells (corresponds to \code{prcomp_result$x}).
#' @param coldata A \code{data.frame} object containing the experimental
#' covariates for which the correlatin is to be determined.
#' @param pcs A numeric vector, containing the corresponding number(s) of PC(s)
#' for which to calculate the correlations. Default: \code{c(1:4)}.
#'
#' @return A \code{data.frame} object with \emph{p} values for each covariate
#' and for each principal component
#'
#' @examples \dontrun{
#' # calculate and plot a bar chart of p-values
#' correlate_PCs(pca_res$cells, data.frame(si[, -"MA_order", with = FALSE]), pcs = 1:4) %>%
#'   melt() %>%
#'   ggplot(., aes(x = Var2, y = -10 * log10(value))) +
#'   geom_bar(stat = "identity", position = position_dodge()) +
#'   facet_wrap(~Var1) +
#'   xlab("")
#' }
#' @export
correlate_PCs <- function(pca_res_samples, coldata, pcs = 1:4) {
  # split the analysis for continuous and categorial
  coldataTypes <- sapply(coldata, class)
  # extract the scores from the pc object
  x <- pca_res_samples

  # make sure pcs don't specify more than what we have
  n.pcs <- min(ncol(x), length(pcs))

  # do it until 1:4 PCs
  res <- matrix(NA, nrow = n.pcs, ncol = ncol(coldata))

  colnames(res) <- colnames(coldata)
  rownames(res) <- paste0("PC", pcs)

  for (i in 1:ncol(res)) {
    # for each covariate...
    for (j in pcs) {
      if (j <= ncol(x)) {
        if (!any(coldataTypes[[i]] %in% c("factor", "character"))) {
          # use correlation for continuous values
          #   print(paste("Correlating PC", pcs[j], "and", names(coldata)[i] ))
          res[j, i] <- cor.test(x[, j], coldata[, i], method = "spearman")$p.value
        } else {
          # had to modify this around to deal with ordered factors etc.
          if (is.null(levels(coldata[, i]))) {
            groups <- factor(coldata[, i])
          } else {
            groups <- coldata[, i]
          }

          if (length(levels(groups)) > 1) {
            res[j, i] <- kruskal.test(x[, j], groups)$p.value
          } else {
            warning(
              paste(
                "Whatever is in",
                colnames(coldata)[i],
                "can neither be used for a Kruskal Test nor a correlation.",
                "Check the content;",
                "it should be either continuous or factors or characters."
              )
            )
          }
        }
      } else {
        warning(
          paste("PC", j, "is not part of the PCA results and will be ignored.")
        )
      }
    }
  }
  return(res)
}

################################################################################
## PCA Plots ###################################################################
################################################################################

#' Plotting PCA plots from DESeq objects
#'
#' @description This function is meant to replace the pcaExplorer::plotPCA function.
#' It borrows heavily from numerous principles established in the \code{scater}
#' package by Aaron Lun.
#'
#' @param object DESeqTransform object; or any object that has its metadata stored
#' in \code{colData(object)}
#' @param pca_results the matrix returned by \code{return_pca_coordinates}
#' @param pc_x the accessor for the coordinates for the x axis of the biplot;
#' default: "PC1"
#' @param pc_y the accessor coordinates for the y axis of the biplot; default:
#' "PC2"
#' @param color_by can be a gene name (part of \code{rownames(object}) or a
#' factor defined in \code{colData(object)}
#' @param shape_by default: \code{NULL}
#' @param size_by default: \code{NULL}
#' @param circle_by Used for encircling points of interest; default: \code{NULL}
#' @param which_quantile only necessary of \code{circle_by} is not NULL. Defines
#' how much the circles should be stretched. The default setting (1) will make sure
#' that all samples per factor are included in the respective circle. If this
#' distorts the circles because of a couple of outlier samples, it may make sense
#' to lower the threshold, e.g. to \code{which_quantile = .99}, which will make
#' the circling ignore crass outliers.
#' @param ... More parameters for \code{plot_pca.df}
#'
#' @return ggplot2 object; either for printing as is or for modifying it even
#' more following the ggplot2 syntax
#'
#' @examples \dontrun{
#' pca_mat <- prep_PCA(assay(myrlt),
#'   n_VarGenes = 500,
#'   input$pca_nrgenes, scale_features = TRUE
#' )
#' pcres <- prcomp(pca_mat)
#' pccoords <- return_PCA_coords(pcres, cell_names = colnames(myrlt), n_PCs = 20)
#'
#' df_to_plot <- data.frame(
#'   x = pccoords$cells[, "PC1"],
#'   y = pccoords$cells[, "PC2"],
#'   row.names = colnames(myrlt)
#' )
#' df_to_plot <- merge(df_to_plot, data.frame(colData(myrlt)), by = "row.names")
#' }
plot_pca_deseq <- function(object, pca_results, pc_x = "PC1", pc_y = "PC2",
                           color_by = NULL, shape_by = NULL, size_by = NULL,
                           circle_by = NULL, which_quantile = 1, ...) {
  # get the data.frame for plotting --------------------------------------------
  ## data.frame with coordinates
  df_to_plot <- data.frame(
    x = pca_results[, pc_x],
    y = pca_results[, pc_y],
    row.names = colnames(object)
  )

  names(df_to_plot) <- c(pc_x, pc_y)
  ## add meta data
  df_to_plot <- merge(df_to_plot, data.frame(colData(object)), by = "row.names")

  # if color/shape aren't already present, add the respective gene values
  if (!is.null(color_by)) {
    df_to_plot <- fx.add_gene_expression(color_by, df_to_plot, object)
  }
  if (!is.null(size_by)) {
    df_to_plot <- fx.add_gene_expression(size_by, df_to_plot, object)
  }

  # get percent variation ----------------------------------------------------
  percentVar <- fx.get_percent_var(pca_results, pcs_select = c(pc_x, pc_y))

  # make the plot ------------------------------------------------------------
  out_plot <- plot_pca.df(df_to_plot,
    color_by = color_by, shape_by = shape_by,
    size_by = size_by, circle_by = circle_by,
    which_quantile = which_quantile,
    percentVar = percentVar, theme_size = 14, legend = "auto",
    alpha = .8, ...
  )

  return(out_plot)
}

#' Wrapper function for plotting PCA plots
#'
#' @param df_to_plot \code{data.frame} for plotting xy plot; should have the
#' values for X and Y in the \strong{first two} columns
#' @param color_by specify the column whose entries will be used to assign either
#' discrete or continous color schemes. Alternatively, you can supply a
#'  \code{list} with TRUE/FALSE logicals and a label stored in \code{list$title}
#' and values stored in \code{list$result} (\code{result} should be as long as
#' \code{dim(object)[2]}.
#' @param shape_by max. 10 levels are accepted for \code{shape_by}; ideally
#'  something of \code{colData}. If you want to specify the (same) shape for
#'  \emph{all} points, use an integer (default shape is 21).
#' @param size_by the factor to assign the size of the points to;
#' default: \code{NULL}. If you want to specify the (same) size for \emph{all}
#' points, use one integer (the default size is 2).
#' @param circle_by a factor (e.g. "condition") that will be used to draw circles
#' around the points belonging to the same factor instance.
#' @param which_quantile only necessary of \code{circle_by} is not NULL. Defines
#' how much the circles should be stretched. The default setting (1) will make sure
#' that all samples per factor are included in the respective circle. If this
#' distorts the circles because of a couple of outlier samples, it may make sense
#' to lower the threshold, e.g. to \code{which_quantile = .99}, which will make
#' the circling ignore crass outliers.
#' @param remove_rug Boolean; default: FALSE.
#' @param percentVar vector of variance values (length should correspond to \code{ncol(df_to_plot)})
#' @param theme_size font size for the plot
#' @param legend choose wether there should be no legend ("none"), a legend for every factor defined by
#' \code{color_by}, \code{shape_by}, \code{size_by} ("all") or whether the legend will only show up for
#' those factors that have more than 1 level ("auto", default setting).
#' @param alpha define the opacity of the points, default: .65
#' @param set_colors Set to \code{FALSE} if you want to add your own color scheme
#' to the resulting plot. The default behavior is to try to automatically assign
#' an optimized coloring scheme.
#' @param set_fill_colors Set to \code{FALSE} if you want to add your own color
#' scheme to the resulting plot \emph{for the background circles}. The default
#' behavior is to try to automatically assign an optimized coloring scheme.
#'
#' @seealso \code{\link{plot_reduced_dim.object}}, which is a wrapper around
#' \code{\link{plot_reduced_dim.df}} and \code{\link{fx.get_reducedDim_df}}
#' @return \code{ggplot2} object (a plot)
#' @examples
#' \dontrun{
#' pca_mat <- prep_PCA(assay(myrlt),
#'   n_VarGenes = 500,
#'   input$pca_nrgenes, scale_features = TRUE
#' )
#' pcres <- prcomp(pca_mat)
#' pccoords <- return_PCA_coords(pcres, cell_names = colnames(myrlt), n_PCs = 20)
#'
#' df_to_plot <- data.frame(
#'   x = pccoords$cells[, "PC1"],
#'   y = pccoords$cells[, "PC2"],
#'   row.names = colnames(myrlt)
#' )
#' df_to_plot <- merge(df_to_plot, data.frame(colData(myrlt)), by = "row.names")
#' }
plot_pca.df <- function(df_to_plot,
                        color_by = NULL, shape_by = NULL,
                        size_by = NULL, circle_by = NULL,
                        which_quantile = 1,
                        remove_rug = FALSE,
                        percentVar = NULL,
                        theme_size = 14, legend = "auto", alpha = .75,
                        set_colors = TRUE, set_fill_colors = TRUE) {
  ## set the defaults for the plot
  if (is.numeric(shape_by)) {
    shape_val <- shape_by
    shape_by <- NULL
  } else {
    shape_val <- 19
  }

  if (is.numeric(size_by)) {
    size_val <- size_by
    size_by <- NULL
  } else {
    size_val <- 3
  }

  ggplot2::update_geom_defaults("point", list(shape = shape_val, size = size_val))

  ## check legend argument
  legend <- match.arg(legend, c("auto", "none", "all"), several.ok = FALSE)

  ## this expects the x and y coordinates to be in the first two columns with name "PC"
  # comps <- colnames(df_to_plot)[1:2]
  comps <- grep("^PC", names(df_to_plot), value = TRUE)[1:2]

  if (is.null(percentVar)) {
    x_lab <- comps[1]
    y_lab <- comps[2]
  } else {
    x_lab <- paste0(
      comps[1], ": ", round(percentVar[1] * 100),
      "% variance"
    )
    y_lab <- paste0(
      comps[2], ": ", round(percentVar[2] * 100),
      "% variance"
    )
  }

  ## generate base plot, apply colour_by, shape_by and size_by variables
  ## (NULL will simply be ignored by ggplot)
  if (!is.null(color_by) && !(color_by %in% names(df_to_plot))) {
    color_by <- NULL
  }
  if (!is.null(shape_by) && !(shape_by %in% names(df_to_plot))) {
    shape_by <- NULL
  }
  if (!is.null(size_by) && !(size_by %in% names(df_to_plot))) {
    size_by <- NULL
  }

  if (is.numeric(df_to_plot[, shape_by])) {
    df_to_plot[, shape_by] <- factor(df_to_plot[, shape_by])
  }

  plot_out <- ggplot(
    df_to_plot,
    aes_string(
      x = comps[1],
      y = comps[2],
      colour = fx.parse_column_names(color_by),
      shape = fx.parse_column_names(shape_by),
      size = fx.parse_column_names(size_by)
    ),
    aes(key = row.names(df_to_plot))
  ) +
    xlab(x_lab) +
    ylab(y_lab)

  if (!remove_rug) {
    plot_out <- plot_out + geom_rug(colour = "gray20", alpha = 0.5)
  }

  plot_out <- plot_out + theme_bw(theme_size)

  ## add background circles
  if (!is.null(circle_by)) {
    c.dt <- fx.get_circle_dt(df_to_plot, fill_by = circle_by, which_quant = which_quantile)
    plot_out <- plot_out +
      ggalt::geom_encircle(
        data = c.dt,
        aes_string(
          fill = fx.parse_column_names(circle_by),
          color = NULL, shape = NULL, size = NULL
        ),
        s_shape = 0.9, expand = 0.07, alpha = .1
      )
  }

  ##  add points
  plot_out <- plot_out + geom_point(alpha = alpha)

  ## try to assign a good color scheme based on the features of color_by
  if (!is.null(color_by) && set_colors) {
    plot_out <- fx.resolve_plot_colours(plot_out, df_to_plot[, color_by],
      colour_by_name = color_by,
      fill = FALSE
    )
  }

  if (!is.null(circle_by) && set_fill_colors) {
    plot_out <- fx.resolve_plot_colours(plot_out, df_to_plot[, circle_by],
      colour_by_name = circle_by,
      fill = TRUE
    )
  }

  ## add sensible shapes
  if (!is.null(shape_by)) {
    plot_out <- fx.resolve_plot_shapes(plot_out, df_to_plot[, shape_by])
  }

  ## Define plotting theme
  plot_out <- plot_out + theme_bw(theme_size)

  ## remove legend if so desired
  if (legend == "none") {
    plot_out <- plot_out + theme(legend.position = "none")
  }

  ## Return plot
  return(plot_out)
}

#' Extract variation captured by PCs
#' @param rd reduced dimension result
#' @param pcs_select specify the principal components that are going to be used
#'
#' @return named vector of variance measurements
fx.get_percent_var <- function(rd, pcs_select = NULL) {
  ## Let's get the variation out of that
  var_vals <- attr(rd, "percentVar")
  if (is.null(var_vals)) {
    return(var_vals)
  } else {
    names(var_vals) <- colnames(rd)

    ## Let's filter it if specific PCs are given
    if (!is.null(pcs_select)) {
      var_vals <- var_vals[pcs_select]
    }
    return(var_vals)
  }
}

#' Make a data.table for encircling points in a reduced dimension plot
#'
#' @description For the encircling, we will use the values of \code{fill_by}
#' to determine the coordinates for the circle. Sometimes, it is useful to remove
#' outlier coordinates before drawing the circles. This is what this function does.
#' If \code{which_quant = 1}, the output will be almost exactly like in.df,
#' except that it will be a data.table.
#'
#' @param in.df data.frame that is used for making the original plot. Should
#' have the x and y coordinates in the first two columns that start with "PC"
#' @param fill_by factor that will be used to determine the background color; if
#' NULL (default), only the outline will be drawn.
#' @param which_quant define the fraction of individual points that should be
#' traced by the circle. Default is to capture all that belong to one instance
#' of \code{fill_by}, i.e. \code{which_quant = 1}. It might be useful to set it
#' to less than 1 (e.g., 0.99) if you don't want the circle to be blown up by
#' individual outliers.
#'
#' @return data.table, most likely a subset of the in.df, depending on the setting
#' of \code{which_quant}. This data.table can be used with ggalt::geom_encircle.
#'
#' @seealso \code{\link{plot_reduced_dim.df}}
#' @examples \dontrun{
#' c.dt <- fx.get_circle_dt(df_to_plot, fill_by = circle_by, which_quant = which_quantile)
#' plot_out <- plot_out + ggalt::geom_encircle(
#'   data = c.dt,
#'   aes_string(fill = fx.parse_column_names(circle_by)),
#'   s_shape = 0.9, expand = 0.07, alpha = .1
#' )
#' }
fx.get_circle_dt <- function(in.df, fill_by = NULL, which_quant = 1) {
  ABCutilities::check_columns(
    fill_by, in.df,
    "data.frame for plotting (df_to_plot)",
    "fx.get_circle_dt"
  )

  comps <- grep("^PC", names(in.df), value = TRUE)[1:2]

  if (!(all(is.character(in.df[, fill_by]))) || is.null(levels(in.df[, fill_by]))) {
    warning(paste("The parameter you chose for the encircling is going to be turned into a factor."))
    in.df[, fill_by] <- factor(in.df[, fill_by])
  }

  pdt <- as.data.table(in.df)
  pdt[, x.high := quantile(get(comps[1]), probs = which_quant), by = fill_by]
  pdt[, y.high := quantile(get(comps[2]), probs = which_quant), by = fill_by]
  pdt[, x.low := quantile(get(comps[1]), probs = 1 - which_quant), by = fill_by]
  pdt[, y.low := quantile(get(comps[2]), probs = 1 - which_quant), by = fill_by]

  out.dt <- pdt[
    get(comps[1]) >= x.low & get(comps[1]) <= x.high &
      get(comps[2]) >= y.low & get(comps[2]) <= y.high,
    colnames(in.df),
    with = FALSE
  ]

  out.df <- as.data.frame(out.dt)

  return(out.df)
}

#' Parse column names
#'
#' @details There are numerous ways to access the values
#' stored in the columns of a \code{data.frame}. The one
#' using a dollar sign (e.g. \code{df$this_column}) is
#' vulnerable to special characters, such as "-", spaces etc.
#' Since the default method of \code{ggplot2} is to use the
#' dollar sign approach, we need to ensure that special
#' characters, which are often part of gene names, are taken
#' care of. This function simply wraps those names into
#' single quotes.
#'
#' @param cn column name
#' @return column name enclosed in single quotes
fx.parse_column_names <- function(cn) {
  if (is.null(cn)) {
    cnp <- cn
  } else {
    cnp <- paste0("`", cn, "`")
  }

  return(cnp)
}

#' Get nice plotting colour schemes for very general colour variables
#' @param plot_out ggplot2 object
#' @param colour_by vector of values that determine the coloring of \code{plot_out}
#' @param colour_by_name string indicating the title/name for \code{colour_by},
#' e.g. the name of a gene
#' @param fill Boolean, default: \code{FALSE}
#' @return \code{ggplot2} object with adjusted coloring scheme
fx.resolve_plot_colours <- function(plot_out, colour_by, colour_by_name,
                                    fill = FALSE) {
  ## if the colour_by object is NULL, return the plot_out object unchanged
  if (is.null(colour_by)) {
    return(plot_out)
  }
  ## Otherwise, set a sensible colour scheme and return the plot_out object
  leg_title <- colour_by_name

  if (fill) { ## routine for fill
    if (is.numeric(colour_by)) {
      plot_out <- plot_out +
        viridis::scale_fill_viridis(name = leg_title)
    } else {
      nlevs_colour_by <- nlevels(as.factor(colour_by))
      if (nlevs_colour_by <= 14) {
        plot_out <- plot_out + scale_fill_manual(
          values = fx.get_palette_ABC("paired_pal"),
          name = leg_title
        )
      } else {
        if (nlevs_colour_by > 14 && nlevs_colour_by <= 20) {
          plot_out <- plot_out + scale_fill_manual(
            values = fx.get_palette_ABC("tableau20"),
            name = leg_title
          )
        } else {
          plot_out <- plot_out +
            viridis::scale_fill_viridis(
              name = leg_title, discrete = TRUE
            )
        }
      }
    }
  } else { ## routine for color
    if (is.numeric(colour_by)) {
      plot_out <- plot_out +
        viridis::scale_color_viridis(name = leg_title)
    } else {
      nlevs_colour_by <- nlevels(as.factor(colour_by))
      if (nlevs_colour_by <= 14) {
        plot_out <- plot_out + scale_colour_manual(
          values = fx.get_palette_ABC("paired_pal"),
          name = leg_title
        )
      } else {
        if (nlevs_colour_by > 14 && nlevs_colour_by <= 20) {
          plot_out <- plot_out + scale_colour_manual(
            values = fx.get_palette_ABC("tableau20"),
            name = leg_title
          )
        } else {
          plot_out <- plot_out +
            viridis::scale_color_viridis(
              name = leg_title, discrete = TRUE
            )
        }
      }
    }
  }
  plot_out
}


#' Color palettes
#'
#' @details Based on \code{scater}'s defaults, but with significant changes to the
#' standard colors that were bing used
fx.get_palette_ABC <- function(palette_name) {
  switch(palette_name,
    # tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    #               "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    #               "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    #               "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"),
    # tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
    #                     "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
    #                     "#CDCC5D", "#6DCCDA"),
    tableau10medium = c(
      "#34B20D", "#FFAE18", "#be2f00", "#73FFC3", "#0D14B2",
      "#8EB20E", "#FF81DE", "#FFFF00", "#0CCC9C", "#656BB2"
    ),
    # from RColorBrewer::brewer.pal(12,  "Paired")
    paired_pal = c(
      "#A6CEE3", "limegreen", "grey30", "grey80", "#1F78B4",
      "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
      "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
    ),
    tableau20 = c(
      "#34B20D", "#FFAE18", "#be2f00", "#73FFC3", "#0D14B2",
      "#8EB20E", "#FF81DE", "#FFFF00", "#0CCC9C", "#656BB2", # tableau10medium
      c(
        "#FF10FC", "#3E68FF", "#8B9440", "#7F85C3", "#FF85C3",
        "#F52306", "#FFD2BD", "#25FFED", "black", "#FFEC83"
      )
    ),
    colorblind10 = c(
      "#006BA4", "#FF800E", "#ABABAB", "#595959",
      "#5F9ED1", "#C85200", "#898989", "#A2C8EC",
      "#FFBC79", "#CFCFCF"
    ),
    trafficlight = c(
      "#B10318", "#DBA13A", "#309343", "#D82526",
      "#FFC156", "#69B764", "#F26C64", "#FFDD71",
      "#9FCD99"
    ),
    purplegray12 = c(
      "#7B66D2", "#A699E8", "#DC5FBD", "#FFC0DA",
      "#5F5A41", "#B4B19B", "#995688", "#D898BA",
      "#AB6AD5", "#D098EE", "#8B7C6E", "#DBD4C5"
    ),
    bluered12 = c(
      "#2C69B0", "#B5C8E2", "#F02720", "#FFB6B0", "#AC613C",
      "#E9C39B", "#6BA3D6", "#B5DFFD", "#AC8763", "#DDC9B4",
      "#BD0A36", "#F4737A"
    ),
    greenorange12 = c(
      "#32A251", "#ACD98D", "#FF7F0F", "#FFB977",
      "#3CB7CC", "#98D9E4", "#B85A0D", "#FFD94A",
      "#39737C", "#86B4A9", "#82853B", "#CCC94D"
    ),
    cyclic = c(
      "#1F83B4", "#1696AC", "#18A188", "#29A03C", "#54A338",
      "#82A93F", "#ADB828", "#D8BD35", "#FFBD4C", "#FFB022",
      "#FF9C0E", "#FF810E", "#E75727", "#D23E4E", "#C94D8C",
      "#C04AA7", "#B446B3", "#9658B1", "#8061B4", "#6F63BB"
    )
  )
}

#' Add gene expression to the data.frame for PCA plotting
#' @param check_factor gene name that will be added
#' @param plot_df data.frame for ggplotting
#' @param ge_object object with gene expression data stored in \code{assay}
#'
#' @return data.frame with additional column for \code{check_factor} values.
fx.add_gene_expression <- function(check_factor, plot_df, ge_object) {
  if (!is.null(check_factor) && !(check_factor %in% names(plot_df))) {
    if (check_factor %in% row.names(ge_object)) {
      out_plot <- merge(plot_df, data.frame(t(assay(ge_object[check_factor, ]))),
        by.x = "Row.names", by.y = "row.names"
      )
    } else {
      warning(paste("Didn't find the factor", check_factor, ". Check it's present in either
                    colData(object) or row.names(object)."))
    }
  } else {
    out_plot <- plot_df
  }

  return(out_plot)
}


#' @details The ggplot2 default setting will complain about having more than 6
#' different factors assigned to shape. This function assigns 12 shapes, trying
#' to find a good order that will still allow for differentiations.
#'
#' @param plot_in ggplot2 object where shape is one of the aesthetics
#' @param shape_by the factors to which shape has been mapped
#'
#' @return ggplot2 object with the shape scale set manually
fx.resolve_plot_shapes <- function(plot_in, shape_by) {
  shps <- c(16, 0, 4, 17, 2, 15, 1, 9, 8, 12, 18, 3)
  fill_na <- length(unique((shape_by))) - length(shps)
  if (fill_na > 0) {
    shps <- c(shps, rep(NA, fill_na))
  }
  plot_out <- plot_in + scale_shape_manual(values = shps)

  return(plot_out)
}


## from scBrowser #############################################################

#' Parse the faceting information
#'
#' @param row_facet name (or NULL) of the factor used for assigning the row facet
#' @param col_facet name (or NULL) of the factor used for assigning the column facet
#' @param wrap Boolean indicating whether facet_grid or facet_wrap should be used.
#' Default: FALSE will return entries fit for facet_grid.
#'
#' @return list with a vector of labels that will be passed to the plotting
#' function (metainfo) and the individual row and column faceting factors
#'
#' @examples \dontrun{
#' fc_info <- choose_faceting(row_facet = input$facet_by_row, col_facet = input$facet_by_col)
#' P <- scABC:::plot_reduced_dim.object(sce,
#'   X = pca$cells[, "PC1"], Y = pca$cells[, "PC2"],
#'   exprs_values = input$expr,
#'   add_cell_info = fc_info$metainfo
#' )
#' gg <- P + facet_grid(as.formula(paste(fc_info$fac_r, "~", fc_info$fac_c)))
#' }
choose_faceting <- function(row_facet, col_facet, wrap = FALSE) {
  if ((row_facet == "none" || is.null(row_facet)) && (col_facet == "none" || is.null(col_facet))) {
    mi <- NULL
    return(list(metainfo = mi))
  } else {
    # decide on row facet
    if (row_facet == "none") {
      fac_r <- NULL
      if (wrap == TRUE) {
        fac_r_out <- NULL
      } else {
        fac_r_out <- "."
      }
    } else {
      fac_r <- row_facet
      fac_r_out <- row_facet
    } # decide on col facet
    if (col_facet == "none") {
      fac_c <- NULL
      if (wrap == TRUE) {
        fac_c_out <- NULL
      } else {
        fac_c_out <- "."
      }
    } else {
      fac_c <- col_facet
      fac_c_out <- col_facet
    }
    mi <- c(fac_r, fac_c)
    return(list(metainfo = mi, fac_r = fac_r_out, fac_c = fac_c_out))
  }
}

#' Select factors of interest based on some customized criteria
#'
#' @description Some of the plotting factors need to be checked for certain
#' features, such as the numbers of levels. For example, you do not want to use
#' a variable with more than 6 levels for assigning the point sizes of a certain
#' plot.
#'
#' @param in.df data.frame variant where the columns represent the different
#' variables/factors that need to be checked, e.g. \code{colData(sce)}.
#' @param max_nlevel integer specifying the maximum number of levels that are
#' allowed. Default: Inf
#' @param factor_class specify the class(es) of the factor that you either want to
#' include or exclude, depending on the setting of \code{keep_class}. Default
#' (NULL) will not check the class.
#' @param keep_class will only come into effect when \code{factor_class} is not
#' NULL; default setting (TRUE) will \strong{retain} all factors that correspond
#' to the class(es) defined by \code{factor_class}. Setting this to FALSE would
#' \strong{remove} those columns.
#'
#' @return The original data.frame without the columns that did not meet the
#' \code{max_nlevel} and/or \code{factor_class} filtering.
#'
#' @examples \dontrun{
#' names(select_factors(colData(sce), max_nlevel = 6))
#' }
select_factors <- function(in.df,
                           max_nlevel = Inf,
                           factor_class = NULL,
                           keep_class = TRUE) {
  # check the levels
  tmp.df <- in.df[, unlist(lapply(in.df, function(x) nlevels(factor(x)) <= max_nlevel))]

  # check factor type if wanted
  if (!is.null(factor_class)) {
    keep_col <- unlist(lapply(tmp.df, function(x) class(x) %in% factor_class))
    if (!keep_class) {
      keep_col <- !keep_col
    }
    tmp.df <- tmp.df[keep_col]
  }
  if (dim(tmp.df)[2] == 0) {
    warning("Select_factors removed all columns from the input data.frame.")
  }
  return(tmp.df)
}
