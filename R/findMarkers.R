# from scran
library(limma)
library(DESeq2)

#attach(loadNamespace("scran"), name = "scran_all")

.subset_to_index <- function(subset, x, byrow=TRUE) {
  if (byrow) {
    dummy <- seq_len(nrow(x))
    names(dummy) <- rownames(x)
  } else {
    dummy <- seq_len(ncol(x))
    names(dummy) <- colnames(x)
  }

  if (!is.null(subset)) {
    dummy <- dummy[subset]
  }
  out <- unname(dummy)
  if (any(is.na(out))) {
    stop("'subset' indices out of range of 'x'")
  }
  return(out)
}

.perform_eb_shrinkage <- function (sigma2, covariate, design)
{
  df.residual <- rep(nrow(design) - ncol(design), length(covariate))
  eb.out <- limma::squeezeVar(sigma2, df = df.residual, robust = TRUE,
                       covariate = covariate)
  df.total <- df.residual + eb.out$df.prior
  df.pooled <- sum(df.residual, na.rm = TRUE)
  df.total <- pmin(df.total, df.pooled)
  return(list(var.post = eb.out$var.post, df.total = df.total))
}

.ranksafe_qr <- function(design, tol=1e-7)
  # Rank-checking QR decomposition of a design matrix. Throws an
  # error if the design matrix is not of full rank, which simplifies
  # downstream processes as full rank can always be assumed.
{
  out <- qr(design, LAPACK=TRUE)
  d <- diag(out$qr)
  if (!all(abs(d) > tol)) {
    stop("design matrix is not of full rank")
  }
  return(out)
}

findMarkers2 <- function (x, clusters, design = NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), min.mean = 0.1, subset.row = NULL)
  {
    clusters <- as.factor(clusters)
    full.design <- model.matrix(~0 + clusters)
    colnames(full.design) <- clust.vals <- levels(clusters)
    if (!is.null(design)) {
      out <- qr.solve(design, cbind(rep(1, nrow(design))))
      to.drop <- abs(out) > 1e-08
      if (any(to.drop)) {
        design <- design[, -which(to.drop)[1], drop = FALSE]
      }
      full.design <- cbind(full.design, design)
    }
    pval.type <- match.arg(pval.type)
    direction <- match.arg(direction)
#    print(direction)
    subset.row <- .subset_to_index(subset.row, x, byrow = TRUE)
    QR <- .ranksafe_qr(full.design)
    stats <- .Call(scran:::cxx_fit_linear_model, QR$qr, QR$qraux,x, subset.row - 1L, TRUE)
    coefficients <- stats[[1]][order(QR$pivot), , drop = FALSE]
    means <- stats[[2]]
    sigma2 <- stats[[3]]
    higher <- means >= min.mean
    if (is.null(min.mean) || all(higher) || !any(higher)) {
      eb.out <- .perform_eb_shrinkage(sigma2, covariate = means,
                                      design = full.design)
    } else {
      high.out <- .perform_eb_shrinkage(sigma2[higher],
                                        covariate = means[higher], design = full.design)
      low.out <- .perform_eb_shrinkage(sigma2[!higher],
                                       covariate = means[!higher], design = full.design)
      eb.out <- mapply(low.out, high.out, FUN = function(low,
                                                         high) {
        val <- numeric(length(means))
        val[higher] <- high
        val[!higher] <- low
        val
      }, SIMPLIFY = FALSE)
    }
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)
    output <- vector("list", length(clust.vals))
    names(output) <- clust.vals
    for (h in seq_along(clust.vals)) {
      host <- clust.vals[h]
      not.h <- seq_along(clust.vals)[-h]
      targets <- clust.vals[not.h]
      con <- matrix(0, ncol(full.design), length(clust.vals))
      diag(con) <- -1
      con[h, ] <- 1
      con <- con[, not.h, drop = FALSE]
      lfit2 <- contrasts.fit(lfit, con)
      ngenes <- length(subset.row)
      ncon <- length(not.h)
      all.lfc <- all.p <- matrix(0, ngenes, ncon)
      ref.coef <- coefficients[h, ]
      for (con in seq_len(ncon)) {
        cur.lfc <- ref.coef - coefficients[not.h[con],
                                           ]
        all.lfc[, con] <- cur.lfc
        cur.t <- cur.lfc/(lfit2$stdev.unscaled[con] *
                            sqrt(eb.out$var.post))
        cur.p <- 2 * pt(-abs(cur.t), df = eb.out$df.total)
        if (direction == "up") {
          cur.p <- ifelse(cur.lfc > 0, cur.p/2, 1 - cur.p/2)
        }
        else if (direction == "down") {
          cur.p <- ifelse(cur.lfc < 0, cur.p/2, 1 - cur.p/2)
        }
        all.p[, con] <- cur.p
      }
      colnames(all.lfc) <- paste0("logFC.", targets)
      if (pval.type == "any") {
        gene.id <- rep(seq_len(ngenes), ncon)
        o <- order(gene.id, all.p)
        penalty <- rep(ncon/seq_len(ncon), ngenes)
        com.p <- matrix(all.p[o] * penalty, ngenes, ncon,
                        byrow = TRUE)
        smallest <- (max.col(-com.p) - 1) * ngenes +
          seq_len(ngenes)
        pval <- com.p[smallest]
      }
      else {
        largest <- (max.col(all.p) - 1) * ngenes + seq_len(ngenes)
        pval <- all.p[largest]
      }
      min.rank <- rep(ngenes, ngenes)
      min.p <- rep(1, ngenes)
      for (con in seq_len(ncon)) {
        cur.p <- all.p[, con]
        cur.rank <- rank(cur.p, ties.method = "first")
        min.rank <- pmin(min.rank, cur.rank)
        min.p <- pmin(min.p, cur.p)
      }
      marker.set <- data.frame(Top = min.rank, Gene = rownames(x)[subset.row],
                               FDR = p.adjust(pval, method = "BH"), all.lfc,
                               stringsAsFactors = FALSE, check.names = FALSE)
      marker.set <- marker.set[order(marker.set$Top, min.p),
                               ]
      rownames(marker.set) <- NULL
      output[[host]] <- marker.set
    }
    return(output)
  }








findMarkers3 <- function (x, clusters, design = NULL, pval.type=c("any", "all"), direction=c("any", "up", "down"), min.mean = 0.1, subset.row = NULL) {
    clusters <- as.factor(clusters)
    full.design <- model.matrix(~0 + clusters)
    colnames(full.design) <- clust.vals <- levels(clusters)
    if (!is.null(design)) {
      out <- qr.solve(design, cbind(rep(1, nrow(design))))
      to.drop <- abs(out) > 1e-08
      if (any(to.drop)) {
        design <- design[, -which(to.drop)[1], drop = FALSE]
      }
      full.design <- cbind(full.design, design)
    }
    pval.type <- match.arg(pval.type)
    direction <- match.arg(direction)
    subset.row <- .subset_to_index(subset.row, x, byrow = TRUE)
    QR <- .ranksafe_qr(full.design)
    stats <- .Call(scran:::cxx_fit_linear_model, QR$qr, QR$qraux,
                   x, subset.row - 1L, TRUE)
    coefficients <- stats[[1]][order(QR$pivot), , drop = FALSE]
    means <- stats[[2]]
    sigma2 <- stats[[3]]
    higher <- means >= min.mean
    if (is.null(min.mean) || all(higher) || !any(higher)) {
      eb.out <- .perform_eb_shrinkage(sigma2, covariate = means,
                                      design = full.design)
    } else {
      high.out <- .perform_eb_shrinkage(sigma2[higher],
                                        covariate = means[higher], design = full.design)
      low.out <- .perform_eb_shrinkage(sigma2[!higher],
                                       covariate = means[!higher], design = full.design)
      eb.out <- mapply(low.out, high.out, FUN = function(low,
                                                         high) {
        val <- numeric(length(means))
        val[higher] <- high
        val[!higher] <- low
        val
      }, SIMPLIFY = FALSE)
    }
    lfit <- lmFit(rbind(seq_len(nrow(full.design))), full.design)
    output <- vector("list", length(clust.vals))
    names(output) <- clust.vals
    for (h in seq_along(clust.vals)) {
      host <- clust.vals[h]
      not.h <- seq_along(clust.vals)[-h]
      targets <- clust.vals[not.h]
      con <- matrix(0, ncol(full.design), length(clust.vals))
      diag(con) <- -1
      con[h, ] <- 1
      con <- con[, not.h, drop = FALSE]
      lfit2 <- contrasts.fit(lfit, con)
      ngenes <- length(subset.row)
      ncon <- length(not.h)
      all.lfc <- all.p <- matrix(0, ngenes, ncon)
      ref.coef <- coefficients[h, ]
      for (con in seq_len(ncon)) {
        cur.lfc <- ref.coef - coefficients[not.h[con],
                                           ]
        all.lfc[, con] <- cur.lfc
        cur.t <- cur.lfc/(lfit2$stdev.unscaled[con] *
                            sqrt(eb.out$var.post))
        cur.p <- 2 * pt(-abs(cur.t), df = eb.out$df.total)
        if (direction == "up") {
          cur.p <- ifelse(cur.lfc > 0, cur.p/2, 1 - cur.p/2)
        }
        else if (direction == "down") {
          cur.p <- ifelse(cur.lfc < 0, cur.p/2, 1 - cur.p/2)
        }
        all.p[, con] <- cur.p
      }
      colnames(all.lfc) <- paste0("logFC.", targets)
      if (pval.type == "any") {
        gene.id <- rep(seq_len(ngenes), ncon)
        o <- order(gene.id, all.p)
        penalty <- rep(ncon/seq_len(ncon), ngenes)
        com.p <- matrix(all.p[o] * penalty, ngenes, ncon,
                        byrow = TRUE)
        smallest <- (max.col(-com.p) - 1) * ngenes +
          seq_len(ngenes)
        pval <- com.p[smallest]
      }
      else {
        largest <- (max.col(all.p) - 1) * ngenes + seq_len(ngenes)
        pval <- all.p[largest]
      }
      min.rank <- rep(ngenes, ngenes)
      min.p <- rep(1, ngenes)
      for (con in seq_len(ncon)) {
        cur.p <- all.p[, con]
        cur.rank <- rank(cur.p, ties.method = "first")
        min.rank <- pmin(min.rank, cur.rank)
        min.p <- pmin(min.p, cur.p)
      }
      marker.set <- data.frame(Top = min.rank, Gene = rownames(x)[subset.row],
                               FDR = p.adjust(pval, method = "BH"), all.lfc,
                               stringsAsFactors = FALSE, check.names = FALSE)
      marker.set <- marker.set[order(marker.set$Top, min.p),]
	  marker.set <- subset(marker.set, Top %in% seq(1:25))
      rownames(marker.set) <- NULL
      output[[host]] <- marker.set
    }
    return(output)
}