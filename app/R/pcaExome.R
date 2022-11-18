pcaplot_exome <- function(x, intgroup = "condition", returnData = FALSE,
                          meta = NULL, title = NULL,
                          pcX = 1, pcY = 2, text_labels = TRUE,
                          point_size = 3,
                          ellipse = TRUE, ellipse.prob = 0.95) {
  pca <- x
  percentVar <- pca$varprop
  #  if (!all(intgroup %in% names(colData(x)))) {
  #      stop("the argument 'intgroup' should specify columns of colData(x)")
  #  }
  intgroup.df <- as.data.frame(meta[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(
    PC1 = pca$eigenvect[, pcX], PC2 = pca$eigenvect[, pcY], group = group,
    intgroup.df, names = pca$sample.id
  )
  colnames(d)[1] <- paste0("PC", pcX)
  colnames(d)[2] <- paste0("PC", pcY)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  d$hjust <- ifelse((sign(d[, paste0("PC", pcX)]) == 1), 0.9, 0.1)

  ## base plot
  g <- ggplot(data = d, aes_string(
    x = paste0("PC", pcX),
    y = paste0("PC", pcY),
    color = "group"
  )) +
    geom_point(size = point_size) +
    xlab(paste0("PC", pcX, ": ", round(percentVar[pcX] * 100, digits = 2), "% variance")) +
    ylab(paste0("PC", pcY, ": ", round(percentVar[pcY] * 100, digits = 2), "% variance"))

  if (ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(d, "group", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x[[paste0("PC", pcX)]], x[[paste0(
        "PC",
        pcY
      )]]))
      mu <- c(mean(x[[paste0("PC", pcX)]]), mean(x[[paste0(
        "PC",
        pcY
      )]]))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2,
        mu,
        FUN = "+"
      ), groups = x$group[1])
    })
    if (nrow(ell) > 0) {
      g <- g + geom_path(data = ell, aes_string(
        x = "X1",
        y = "X2", color = "groups", group = "groups"
      ))
    }
  }
  if (text_labels) {
    g <- g + geom_label_repel(mapping = aes_string(
      label = "names",
      fill = "group"
    ), color = "white", show.legend = TRUE) +
      theme_bw()
  }
  if (!is.null(title)) {
    g <- g + ggtitle(title)
  }
  g
}













read.maf.2 <- function(maf, clinicalData = NULL, removeDuplicatedVariants = TRUE,
                       useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
                       gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = "all",
                       cnTable = NULL, isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE) {
  if (is.data.frame(x = maf)) {
    maf <- data.table::setDT(maf)
  } else {
    message("reading maf..")
    if (as.logical(length(grep(
      pattern = "gz$", x = maf,
      fixed = FALSE
    )))) {
      if (Sys.info()[["sysname"]] == "Windows") {
        maf.gz <- gzfile(description = maf, open = "r")
        suppressWarnings(maf <- data.table::as.data.table(read.csv(
          file = maf.gz,
          header = TRUE, sep = "\t", stringsAsFactors = FALSE,
          comment.char = "#"
        )))
        close(maf.gz)
      } else {
        maf <- suppressWarnings(data.table::fread(
          input = paste(
            "zcat <",
            maf
          ), sep = "\t", stringsAsFactors = FALSE,
          verbose = FALSE, data.table = TRUE, showProgress = TRUE,
          header = TRUE
        ))
      }
    } else {
      suppressWarnings(maf <- data.table::fread(
        input = maf,
        sep = "\t", stringsAsFactors = FALSE, verbose = FALSE,
        data.table = TRUE, showProgress = TRUE, header = TRUE
      ))
    }
  }
  maf <- maftools:::validateMaf(
    maf = maf, isTCGA = isTCGA, rdup = removeDuplicatedVariants,
    chatty = verbose
  )
  if (!useAll) {
    message("--Using only `Somatic` variants from Mutation_Status. Set useAll = TRUE to include everything.")
    if (length(colnames(maf)[colnames(x = maf) %in% "Mutation_Status"]) >
      0) {
      maf <- maf[Mutation_Status %in% "Somatic"]
      if (nrow(maf) == 0) {
        stop("No more Somatic mutations left after filtering for Mutation_Status! Maybe set useAll to TRUE ?")
      }
    } else {
      message("---Oops! Mutation_Status not found. Assuming all variants are Somatic and validated.")
    }
  }
  if (is.null(vc_nonSyn)) {
    vc.nonSilent <- c(
      "Frame_Shift_Del", "Frame_Shift_Ins",
      "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation",
      "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
      "Missense_Mutation"
    )
  } else {
    vc.nonSilent <- vc_nonSyn
  }
  maf.silent <- maf[!Variant_Classification %in% vc.nonSilent]
  if (nrow(maf.silent) > 0) {
    maf.silent.vc <- maf.silent[, .N, .(
      Tumor_Sample_Barcode,
      Variant_Classification
    )]
    maf.silent.vc.cast <- data.table::dcast(
      data = maf.silent.vc,
      formula = Tumor_Sample_Barcode ~ Variant_Classification,
      fill = 0, value.var = "N"
    )
    summary.silent <- data.table::data.table(
      ID = c(
        "Samples",
        colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]
      ),
      N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,
        2:ncol(maf.silent.vc.cast),
        with = FALSE
      ]))
    )
    maf <- maf[Variant_Classification %in% vc.nonSilent]
    if (verbose) {
      message(paste0("silent variants: ", nrow(maf.silent)))
      print(summary.silent)
    }
  }
  if (!is.null(gisticAllLesionsFile)) {
    gisticIp <- readGistic(
      gisticAllLesionsFile = gisticAllLesionsFile,
      gisticAmpGenesFile = gisticAmpGenesFile, gisticDelGenesFile = gisticDelGenesFile,
      isTCGA = isTCGA, gisticScoresFile = gisticScoresFile,
      cnLevel = cnLevel
    )
    gisticIp <- gisticIp@data
    suppressWarnings(gisticIp[, `:=`(id, paste(Hugo_Symbol,
      Tumor_Sample_Barcode,
      sep = ":"
    ))])
    gisticIp <- gisticIp[!duplicated(id)]
    gisticIp[, `:=`(id, NULL)]
    maf <- rbind(maf, gisticIp, fill = TRUE)
  } else if (!is.null(cnTable)) {
    if (verbose) {
      message("Processing copy number data..")
    }
    if (is.data.frame(cnTable)) {
      cnDat <- data.table::setDT(cnTable)
    } else {
      cnDat <- data.table::fread(
        input = cnTable, sep = "\t",
        stringsAsFactors = FALSE, header = TRUE, colClasses = "character"
      )
    }
    colnames(cnDat) <- c(
      "Hugo_Symbol", "Tumor_Sample_Barcode",
      "Variant_Classification"
    )
    if (isTCGA) {
      cnDat[, `:=`(Tumor_Sample_Barcode, substr(
        x = cnDat$Tumor_Sample_Barcode,
        start = 1, stop = 12
      ))]
    }
    cnDat$Variant_Type <- "CNV"
    suppressWarnings(cnDat[, `:=`(id, paste(Hugo_Symbol,
      Tumor_Sample_Barcode,
      sep = ":"
    ))])
    cnDat <- cnDat[!duplicated(id)]
    cnDat[, `:=`(id, NULL)]
    maf <- rbind(maf, cnDat, fill = TRUE)
  }
  maf$Variant_Type <- as.factor(as.character(maf$Variant_Type))
  maf$Variant_Classification <- as.factor(as.character(maf$Variant_Classification))
  maf$Tumor_Sample_Barcode <- as.factor(as.character(maf$Tumor_Sample_Barcode))
  if (verbose) {
    message("Summarizing..")
  }
  message("Done !")
  return(maf)
}
