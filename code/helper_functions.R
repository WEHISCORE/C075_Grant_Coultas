# Helper function to coerce a DataFrame to a data.frame while preserving column
# names as is.
.adf <- function(x) {
  setNames(as.data.frame(x), colnames(x))
}

# Helper function to Combine data from 2 SCEs using gene names.
# NOTE: This assumes more than I'd like about the rowData and doesn't do much
#       checking of these assumptions.
.combine <- function(x, y, rowData_by = c("ENSEMBL", "SYMBOL", "CHR")) {
  if (is.null(rowData_by)) {
    rowData <- dplyr::full_join(
      .adf(rowData(x)) %>%
        tibble::rownames_to_column(var = "gene"),
      .adf(rowData(y)) %>%
        tibble::rownames_to_column(var = "gene")) %>%
      tibble::column_to_rownames("gene") %>%
      DataFrame(., row.names = rownames(.))
  } else {
    rowData <- dplyr::full_join(
      .adf(rowData(x)[, rowData_by, drop = FALSE]),
      .adf(rowData(y)[, rowData_by, drop = FALSE]),
      by = rowData_by) %>%
      DataFrame(row.names = scater::uniquifyFeatureNames(
        .$ENSEMBL,
        .$SYMBOL))
    rownames(x) <- rownames(rowData)[match(rowData(x)$ENSEMBL, rowData$ENSEMBL)]
    rownames(y) <- rownames(rowData)[match(rowData(y)$ENSEMBL, rowData$ENSEMBL)]
  }

  colData <- rbind(colData(x), colData(y))

  counts <- matrix(
    data = 0L,
    nrow = nrow(rowData), ncol = nrow(colData),
    dimnames = list(rownames(rowData), rownames(colData)))
  counts[rownames(x), colnames(x)] <- counts(
    x,
    withDimnames = FALSE)
  counts[rownames(y), colnames(y)] <- counts(
    y,
    withDimnames = FALSE)

  stopifnot(
    identical(
      metadata(x)$scPipe$version,
      metadata(y)$scPipe$version))
  stopifnot(
    identical(
      metadata(x)$scPipe$QC_cols,
      metadata(y)$scPipe$QC_cols))
  stopifnot(
    identical(
      metadata(x)$scPipe$demultiplex_info$status,
      metadata(y)$scPipe$demultiplex_info$status))
  stopifnot(
    identical(
      metadata(x)$scPipe$UMI_dup_info$duplication.number,
      metadata(y)$scPipe$UMI_dup_info$duplication.number))
  stopifnot(identical(metadata(x)$Biomart, metadata(y)$Biomart))
  metadata <- list(
    scPipe = list(
      version = metadata(x)$scPipe$version,
      QC_cols = metadata(x)$scPipe$QC_cols,
      demultiplex_info = data.frame(
        status = metadata(x)$scPipe$demultiplex_info$status,
        count = metadata(x)$scPipe$demultiplex_info$count +
          metadata(y)$scPipe$demultiplex_info$count),
      UMI_dup_info = data.frame(
        duplication.number = metadata(
          x)$scPipe$UMI_dup_info$duplication.number,
        count = metadata(x)$scPipe$UMI_dup_info$count +
          metadata(y)$scPipe$UMI_dup_info$count)),
    Biomart = metadata(x)$Biomart)

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(counts = counts),
    metadata = metadata)

  stopifnot(identical(int_metadata(x), int_metadata(y)))
  int_metadata(sce) <- int_metadata(x)
  int_elementMetadata <- dplyr::full_join(
    x = .adf(int_elementMetadata(x)) %>%
      tibble::add_column(gene = rownames(x)),
    y = .adf(int_elementMetadata(y)) %>%
      tibble::add_column(gene = rownames(y))) %>%
    tibble::column_to_rownames("gene") %>%
    DataFrame()
  int_elementMetadata(sce) <- int_elementMetadata

  stopifnot(validObject(sce))
  sce
}

.cbindSCEs <- function(list_of_sce, rowData_by = 1:6) {
  do.call(
    cbind,
    lapply(list_of_sce, function(sce) {
      # NOTE: Some fudging to combine only the necessary bits of each SCE
      #       (basically, don't include any QC metrics).
      rowData(sce) <- rowData(sce)[, rowData_by]
      sce
    }))
}

# NOTE: My best guess of what this function does https://github.com/MarioniLab/compareSingleCell/blob/543aa28e3ae25fad4ffb1d47c27c8a364966095c/vignettes/embryo_expression.Rmd#L76
.sumCountsAcrossCells <- function(sce, cluster_sample) {
  counts <- counts(sce, withDimnames = FALSE)
  edgeR::sumTechReps(counts, cluster_sample)
}

# A modified version of scater::plotHeatmap() that allows me to pass
# `annotation_colors` down to pheatmap::pheatmap().
.plotHeatmap <- function (object, features, columns = NULL, exprs_values = "logcounts",
                          center = FALSE, zlim = NULL, symmetric = FALSE, color = NULL,
                          colour_columns_by = NULL, order_columns_by = NULL, by_exprs_values = exprs_values,
                          by_show_single = FALSE, show_colnames = FALSE, cluster_cols = is.null(order_columns_by),
                          ...) {
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  heat.vals <- assay(object, exprs_values)[features, , drop = FALSE]
  if (!is.null(columns)) {
    columns <- scater:::.subset2index(columns, object, byrow = FALSE)
    heat.vals <- heat.vals[, columns, drop = FALSE]
  }
  if (center) {
    heat.vals <- heat.vals - DelayedMatrixStats::rowMeans2(DelayedArray(heat.vals))
  }
  if (!is.null(order_columns_by)) {
    ordering <- list()
    for (i in seq_along(order_columns_by)) {
      vals <- retrieveCellInfo(object, order_columns_by[[i]],
                               exprs_values = by_exprs_values)$val
      if (!is.null(columns)) {
        vals <- vals[columns]
      }
      ordering[[i]] <- vals
    }
    heat.vals <- heat.vals[, do.call(order, ordering), drop = FALSE]
    cluster_cols <- FALSE
    colour_columns_by <- c(colour_columns_by, order_columns_by)
  }
  if (is.null(zlim)) {
    zlim <- range(heat.vals)
  }
  if (symmetric) {
    extreme <- max(abs(zlim))
    zlim <- c(-extreme, extreme)
  }
  heat.vals[heat.vals < zlim[1]] <- zlim[1]
  heat.vals[heat.vals > zlim[2]] <- zlim[2]
  if (is.null(color)) {
    color <- eval(formals(pheatmap::pheatmap)$color, envir = environment(pheatmap::pheatmap))
  }
  color.breaks <- seq(zlim[1], zlim[2], length.out = length(color) +
                        1L)
  if (length(colour_columns_by)) {
    column_variables <- column_colorings <- list()
    for (i in seq_along(colour_columns_by)) {
      field <- colour_columns_by[[i]]
      colour_by_out <- retrieveCellInfo(object, field,
                                        exprs_values = by_exprs_values)
      if (is.null(colour_by_out$val)) {
        next
      } else if (is.numeric(colour_by_out$val)) {
        colour_fac <- colour_by_out$val
        col_scale <- viridis(25)
      } else {
        colour_fac <- as.factor(colour_by_out$val)
        nlevs_colour_by <- nlevels(colour_fac)
        # if (nlevs_colour_by <= 10) {
        #   col_scale <- .get_palette("tableau10medium")
        # }
        # else if (nlevs_colour_by > 10 && nlevs_colour_by <=
        #          20) {
        #   col_scale <- scater:::.get_palette("tableau20")
        # }
        # else {
        #   col_scale <- viridis(nlevs_colour_by)
        # }
        # col_scale <- col_scale[seq_len(nlevs_colour_by)]
        # names(col_scale) <- levels(colour_fac)
      }
      col_name <- colour_by_out$name
      if (col_name == "") {
        col_name <- paste0("unnamed", i)
      }
      column_variables[[col_name]] <- colour_fac
      # column_colorings[[col_name]] <- col_scale
    }
    column_variables <- do.call(data.frame, c(column_variables,
                                              list(row.names = colnames(object))))
  } else {
    column_variables <- column_colorings <- NULL
  }
  pheatmap::pheatmap(heat.vals, color = color, breaks = color.breaks,
                     annotation_col = column_variables,
                     # annotation_colors = column_colorings,
                     show_colnames = show_colnames, cluster_cols = cluster_cols,
                     ...)
}

# Plot SingleR scores on a reduced dimension plot from a SCE.
plotScoreReducedDim <- function(results, sce, dimred = "TSNE",
                                max.labels = 20, normalize = TRUE, ncol = 5,
                                ...) {
  scores <- results$scores
  rownames(scores) <- rownames(results)
  m <- rowMaxs(scale(t(scores)))
  to.keep <- head(order(m, decreasing = TRUE), max.labels)
  if (normalize) {
    mmax <- rowMaxs(scores)
    mmin <- rowMins(scores)
    scores <- (scores - mmin) / (mmax - mmin)
    scores <- scores ^ 3
  }
  scores <- scores[, to.keep, drop = FALSE]
  cns <- colnames(scores)
  p <- lapply(cns, function(cn) {
    scater::plotReducedDim(
      sce,
      dimred = dimred,
      colour_by = data.frame(Score = scores[, cn]),
      ...) +
      ggtitle(cn) +
      scale_fill_viridis_c(limits = force(if(normalize) c(0, 1) else NULL)) +
      guides(fill = guide_colourbar(title = "Score"))
  })
  cowplot::plot_grid(plotlist = p, ncol = ncol)
}

# Take a DataFrame with AtomicList columns and return a DataFrame where these
# columns have been flattened by paste-ing together the elements separated by
# `sep`.
flattenDF <- function(x, sep = "; ") {
  endoapply(x, function(xx) {
    if (!is(xx, "AtomicList")) {
      return(xx)
    }
    unstrsplit(as(xx, "CharacterList"), sep = sep)
  })
}
