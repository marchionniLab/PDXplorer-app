# library(RCircos)
# library(readr)
# library(chimeraviz)
importStarfusion2 <- function(filename, genomeVersion, limit) {
  validGenomes <- c("hg19", "hg38")
  if (is.na(match(tolower(genomeVersion), tolower(validGenomes)))) {
    stop("Invalid genome version given")
  }
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }
  report <- tryCatch(
    {
      col_types_starfusion <- readr::cols_only(
        `#FusionName` = col_skip(),
        JunctionReadCount = col_integer(),
        SpanningFragCount = col_integer(),
        SpliceType = col_skip(),
        LeftGene = col_character(),
        LeftBreakpoint = col_character(),
        RightGene = col_character(),
        RightBreakpoint = col_character(),
        LargeAnchorSupport = col_character(),
        LeftBreakDinuc = col_character(),
        LeftBreakEntropy = col_number(),
        RightBreakDinuc = col_character(),
        RightBreakEntropy = col_number()
      )
      if (missing(limit)) {
        readr::read_tsv(file = filename, col_types = col_types_starfusion)
      } else {
        readr::read_tsv(
          file = filename, col_types = col_types_starfusion,
          n_max = limit
        )
      }
    },
    error = function(cond) {
      message("Reading filename caused an error:")
      stop(cond)
    },
    warning = function(cond) {
      message("Reading filename caused a warning:")
      warning(cond)
    }
  )
  id <- NA
  inframe <- NA
  fusionTool <- "starfusion"
  spanningReadsCount <- NA
  splitReadsCount <- NA
  junctionSequence <- NA
  fusionReadsAlignment <- NA
  fusionList <- vector("list", dim(report)[1])
  for (i in 1:dim(report)[1]) {
    fusionToolSpecificData <- list()
    fusionToolSpecificData[["LargeAnchorSupport"]] <- report[[
      i,
      "LargeAnchorSupport"
    ]]
    fusionToolSpecificData[["LeftBreakDinuc"]] <- report[[
      i,
      "LeftBreakDinuc"
    ]]
    fusionToolSpecificData[["LeftBreakEntropy"]] <- report[[
      i,
      "LeftBreakEntropy"
    ]]
    fusionToolSpecificData[["RightBreakDinuc"]] <- report[[
      i,
      "RightBreakDinuc"
    ]]
    fusionToolSpecificData[["RightBreakEntropy"]] <- report[[
      i,
      "RightBreakEntropy"
    ]]
    id <- as.character(i)
    leftBreakPoint <- unlist(strsplit(report[[i, "LeftBreakpoint"]],
      split = ":"
    ))
    rightBreakPoint <- unlist(strsplit(report[[i, "RightBreakpoint"]],
      split = ":"
    ))
    chromosomeA <- leftBreakPoint[1]
    chromosomeB <- rightBreakPoint[1]
    breakpointA <- as.numeric(leftBreakPoint[2])
    breakpointB <- as.numeric(rightBreakPoint[2])
    strandA <- leftBreakPoint[3]
    strandB <- rightBreakPoint[3]
    splitReadsCount <- report[[i, "JunctionReadCount"]]
    spanningReadsCount <- report[[i, "SpanningFragCount"]]
    junctionSequenceA <- Biostrings::DNAString()
    junctionSequenceB <- Biostrings::DNAString()
    geneNames1 <- unlist(strsplit(report[[i, "LeftGene"]],
      split = "\\\\^"
    ))
    geneNames1 <- gsub("\\^.*", "", geneNames1)

    geneNames2 <- unlist(strsplit(report[[i, "RightGene"]],
      split = "\\\\^"
    ))
    geneNames2 <- gsub("\\^.*", "", geneNames2)
    nameA <- geneNames1[1]
    nameB <- geneNames2[1]
    ensemblIdA <- NA_character_
    ensemblIdB <- NA_character_
    geneA <- new(
      Class = "PartnerGene",
      name = nameA,
      ensemblId = ensemblIdA,
      chromosome = chromosomeA,
      breakpoint = breakpointA,
      strand = strandA,
      junctionSequence = junctionSequenceA,
      transcripts = GenomicRanges::GRangesList()
    )
    geneB <- new(
      Class = "PartnerGene",
      name = nameB,
      ensemblId = ensemblIdB,
      chromosome = chromosomeB,
      breakpoint = breakpointB,
      strand = strandB,
      junctionSequence = junctionSequenceB,
      transcripts = GenomicRanges::GRangesList()
    )
    fusionList[[i]] <- new(
      Class = "Fusion",
      id = id, fusionTool = fusionTool,
      genomeVersion = genomeVersion,
      spanningReadsCount = spanningReadsCount,
      splitReadsCount = splitReadsCount,
      fusionReadsAlignment = Gviz::AlignmentsTrack(),
      geneA = geneA,
      geneB = geneB,
      inframe = inframe,
      fusionToolSpecificData = fusionToolSpecificData
    )
  }
  fusionList
}


.scaleListToInterval <- function(theList, newMin, newMax) {
  if (length(theList) <= 1) {
    stop(paste(
      "Invalid list. Using this function with less than two values",
      "makes to sense."
    ))
  }
  (newMax - newMin) * (theList - min(theList)) / (max(theList) - min(theList)) + newMin
}

.fusionsToLinkData <- function(fusionList, minLinkWidth = 1, maxLinkWidt = 10) {
  Chromosome <- vector(mode = "character", length = length(fusionList))
  chromStart <- vector(mode = "numeric", length = length(fusionList))
  chromEnd <- vector(mode = "numeric", length = length(fusionList))

  Chromosome.1 <- vector(mode = "character", length = length(fusionList))
  chromStart.1 <- vector(mode = "numeric", length = length(fusionList))
  chromEnd.1 <- vector(mode = "numeric", length = length(fusionList))

  linkWidth <- vector(mode = "numeric", length = length(fusionList))

  for (i in 1:length(fusionList)) {
    fusion <- fusionList[[i]]

    Chromosome[[i]] <- fusion@geneA@chromosome
    chromStart[[i]] <- fusion@geneA@breakpoint
    chromEnd[[i]] <- fusion@geneA@breakpoint + 1 # This value shouldn't matter

    Chromosome.1[[i]] <- fusion@geneB@chromosome
    chromStart.1[[i]] <- fusion@geneB@breakpoint
    chromEnd.1[[i]] <- fusion@geneB@breakpoint + 1 # This value shouldn't matter

    # Set link width = number of spanning+split reads supporting fusion
    linkWidth[[i]] <- fusion@spanningReadsCount + fusion@splitReadsCount
  }

  # Normalize all link width values to the interval [minLinkWidth, maxLinkWidt]
  linkWidth <- .scaleListToInterval(linkWidth, minLinkWidth, maxLinkWidt)

  data.frame(
    Chromosome,
    chromStart,
    chromEnd,
    Chromosome.1,
    chromStart.1,
    chromEnd.1,
    linkWidth
  )
}

.fusionsToGeneLabelData <- function(fusionList) {
  originalLength <- length(fusionList)
  newLength <- originalLength * 2

  Chromosome <- vector(mode = "character", length = newLength)
  chromStart <- vector(mode = "numeric", length = newLength)
  chromEnd <- vector(mode = "numeric", length = newLength)
  Gene <- vector(mode = "character", length = newLength)

  for (i in 1:length(fusionList)) {
    fusion <- fusionList[[i]]

    # We are building the list of gene names for both partner genes at once, so
    # let's put the gene1 genes in the first half of the list..
    Chromosome[[i]] <- fusion@geneA@chromosome
    chromStart[[i]] <- fusion@geneA@breakpoint
    chromEnd[[i]] <- fusion@geneA@breakpoint + 1 # This value shouldn't matter
    Gene[[i]] <- fusion@geneA@name

    # and put the gene2 genes in the other half.
    Chromosome[[i + originalLength]] <- fusion@geneB@chromosome
    chromStart[[i + originalLength]] <- fusion@geneB@breakpoint
    chromEnd[[i + originalLength]] <- fusion@geneB@breakpoint + 1 # This value shouldn't matter
    Gene[[i + originalLength]] <- fusion@geneB@name
  }

  data.frame(
    Chromosome,
    chromStart,
    chromEnd,
    Gene
  )
}


plotCircle2 <- function(fusionList, cytoband) {
  if (!is.list(fusionList)) {
    stop("Input must be a list of fusion objects")
  }
  if (class(fusionList[[1]]) != "Fusion") {
    stop("items in fusionList must be an object of type Fusion")
  }
  cytobandFile <- cytoband
  cytoband <- utils::read.table(cytobandFile)
  names(cytoband) <- c(
    "Chromosome", "ChromStart", "ChromEnd",
    "Band", "Stain"
  )
  assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)

  cytoband <- RCircos.Sort.Genomic.Data(
    genomic.data = cytoband,
    is.ideo = TRUE
  )
  cyto.info <- cytoband
  chr.exclude <- NULL
  tracks.inside <- 3
  tracks.outside <- 0
  RCircos::RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
  params <- RCircos::RCircos.Get.Plot.Parameters()
  params$text.size <- 1
  params$char.width <- 1000
  params$text.color <- "black"
  RCircos.Reset.Plot.Parameters(params)
  RCircos.Env$RCircos.PlotPar$text.size <- .6
  RCircos::RCircos.Set.Plot.Area(margins = 0)
  RCircos::RCircos.Chromosome.Ideogram.Plot()
  geneLabelData <- .fusionsToGeneLabelData(fusionList)
  geneLabelData <- geneLabelData[grep("Un", geneLabelData$Chromosome, invert = T), ]
  geneLabelData <- geneLabelData[grep("converted", geneLabelData$Chromosome, invert = T), ]
  geneLabelData$Gene <- gsub("\\^.*", "", geneLabelData$Gene)
  name.col <- 4
  side <- "in"
  track.num <- 1
  RCircos::RCircos.Gene.Connector.Plot(geneLabelData, track.num, side)
  track.num <- 2
  RCircos::RCircos.Gene.Name.Plot(
    geneLabelData, name.col, track.num,
    side
  )
  linkData <- .fusionsToLinkData(fusionList)
  linkData <- linkData[grep("Un", linkData$Chromosome, invert = T), ]
  linkData <- linkData[grep("converted", linkData$Chromosome, invert = T), ]
  linkData <- linkData[grep("Un", linkData$Chromosome.1, invert = T), ]
  linkData <- linkData[grep("converted", linkData$Chromosome.1, invert = T), ]

  track.num <- 3
  RCircos::RCircos.Link.Plot(
    link.data = linkData, track.num = track.num,
    by.chromosome = TRUE, start.pos = NULL, genomic.columns = 3,
    is.sorted = FALSE, lineWidth = linkData$linkWidth
  )
  remove("RCircos.Env", envir = .GlobalEnv)
}
