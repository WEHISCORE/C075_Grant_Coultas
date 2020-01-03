library(here)
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.scPipe.SCE.rds"))

library(scater)
library(dplyr)
library(janitor)

as.data.frame(colData(sce)) %>%
  group_by(plate_number, cell_type_descriptor, sample_name) %>%
  count() %>%
  knitr::kable()

library(EnsDb.Mmusculus.v79)
ensembl <- gsub("\\.[0-9]+$", "", rownames(sce))
symb <- mapIds(
  x = EnsDb.Mmusculus.v79,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = ensembl,
  keytype = "GENEID",
  column = "SYMBOL")
rowData(sce)$ENSEMBL <- ensembl
rowData(sce)$SYMBOL <- symb
# Replace the row names of the SCE by the gene symbols (where available).
library(scater)
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$SYMBOL)
# Add chromosome location so we can filter on mitochondrial genes.
location <- mapIds(
  x = EnsDb.Mmusculus.v79,
  # NOTE: Need to remove gene version number prior to lookup.
  keys = rowData(sce)$ENSEMBL,
  keytype = "GENEID",
  column = "SEQNAME")
rowData(sce)$CHR <- location

is_mito <- rowData(sce)$CHR == "MT"
sce <- addPerCellQC(sce, subsets = list(Mt = which(is_mito)))

sce$sample <- paste0(
  sce$plate_number,
  ".",
  stringr::str_extract(sce$cell_type_descriptor, "# [0-9]+"))
sce$group <- ifelse(grepl("Mutant", sce$cell_type_descriptor), "mutant", "control")
plotColData(
  sce,
  y = "sum",
  x = "sample",
  colour_by = "group",
  other_fields = "sample_name") +
  scale_y_log10() +
  facet_grid(~sample_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
plotColData(
  sce,
  y = "detected",
  x = "sample",
  colour_by = "group",
  other_fields = "sample_name") +
  scale_y_log10() +
  facet_grid(~sample_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
plotColData(
  sce,
  y = "altexps_ERCC_percent",
  x = "sample",
  colour_by = "group",
  other_fields = "sample_name") +
  ylim(0, 100) +
  facet_grid(~sample_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
plotColData(
  sce,
  y = "subsets_Mt_percent",
  x = "sample",
  colour_by = "group",
  other_fields = "sample_name") +
  ylim(0, 100) +
  facet_grid(~sample_name) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))

p <- lapply(sort(unique(sce$plate_number)), function(p) {
  z <- sce[, sce$plate_number == p]
  plotPlatePosition(
    z,
    as.character(z$well_position),
    colour_by = "altexps_ERCC_percent",
    point_size = 2,
    point_alpha = 1,
    theme_size = 5,
    shape_by = "group",
    by_show_single = TRUE) +
    ggtitle(p) +
    scale_colour_viridis_c(
      limits = c(0, 60),
      breaks = seq(0, 60, 10)) +
    scale_shape_manual(values = c(control = 16, mutant = 17)) +
    guides(shape = FALSE)
})
multiplot(plotlist = p, cols = 3)

# Check expression of mutant gene (Kat7)
plotExpression(
  sce,
  "Kat7",
  x = "sample",
  exprs_values = "counts",
  colour_by = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5))
table(sce$cell_type_descriptor, counts(sce)["Kat7", ] > 0)

sce <- sce[, sce$sample_name == "single cell"]
sce$sample <- factor(sce$sample)
sce$sample_name <- factor(sce$sample_name)

libsize_drop <- isOutlier(
  metric = sce$sum,
  nmads = 3,
  type = "lower",
  log = TRUE)
feature_drop <- isOutlier(
  metric = sce$detected,
  nmads = 3,
  type = "lower",
  log = TRUE)
spike_drop <- isOutlier(
  metric = sce$altexps_ERCC_percent,
  nmads = 5,
  type = "higher")
sce_pre_QC_outlier_removal <- sce
sce <- sce[, !(libsize_drop | feature_drop | spike_drop)]
data.frame(
  ByLibSize = tapply(
    libsize_drop,
    sce_pre_QC_outlier_removal$plate_number,
    sum),
  ByFeature = tapply(
    feature_drop,
    sce_pre_QC_outlier_removal$plate_number,
    sum),
  BySpike = tapply(
    spike_drop,
    sce_pre_QC_outlier_removal$plate_number,
    sum),
  Remaining = as.vector(unname(table(sce$plate_number)))) %>%
  knitr::kable(
    caption = "Number of samples removed by each QC step and the number of samples remaining.")
cbind(
  "pre-QC" = table(sce_pre_QC_outlier_removal$cell_type_descriptor),
  "post-QC" = table(sce$cell_type_descriptor))

plotHighestExprs(sce, n = 50)
plotHighestExprs(sce[!grepl("^mt|^Rpl|^Rps", rownames(sce))], n = 50)

ave_counts <- calcAverage(sce, use_size_factors = FALSE)
par(mfrow = c(1, 1))
hist(
  x = log10(ave_counts),
  breaks = 100,
  main = "",
  col = "grey",
  xlab = expression(Log[10] ~ "average count"))
# NOTE: I've opted to filter out genes with zero counts across all data sets
#       rather than a per-data set basis. This means that each data set continues
#       to have the same set of genes after this filter is applied.
to_keep <- ave_counts > 0
table(to_keep)
sce <- sce[to_keep, ]

library(scran)
set.seed(1011220)
clusters <- quickCluster(
  sce,
  min.size = 60,
  use.ranks = FALSE,
  BSPARAM = BiocSingular::IrlbaParam())
table(clusters)
sce <- computeSumFactors(
  sce,
  clusters = clusters,
  min.mean = 0.1)
summary(sizeFactors(sce))

xlim <- c(1, max(sce$sum) / 1e3)
ylim <- range(sizeFactors(sce))
par(mfrow = c(3, 3))
lapply(sort(unique(sce$plate_number)), function(p) {
  sce <- sce[, sce$plate_number == p]
  plot(
    x = sce$sum / 1e3,
    y = sizeFactors(sce),
    log = "xy",
    xlab = "Library size (thousands)",
    ylab = "Size factor",
    main = p,
    xlim = xlim,
    ylim = ylim,
    pch = 16)
})

sce <- computeSpikeFactors(sce, spikes = "ERCC", general.use = FALSE)
sce <- logNormCounts(sce)

var.out <- modelGeneVarWithSpikes(sce, "ERCC")
par(mfrow = c(1, 1))
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
fit <- metadata(var.out)
points(fit$mean, fit$var, col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes])
set.seed(1000)
sce <- denoisePCA(sce, technical=fit$trend, BSPARAM=BiocSingular::IrlbaParam())
ncol(reducedDim(sce, "PCA"))

set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50)
plotTSNE(sce, colour_by = "group")
plotTSNE(sce, colour_by = "cell_type_descriptor")
plotTSNE(sce, colour_by = "plate_number")
plotTSNE(sce, colour_by = "subsets_Mt_percent")

# TODO: Figure out if its MT, X, and/or Y that's causing the clustering.
# NOTE: Remove MT, X, and Y and re-do t-SNE.
sce <- sce[which(!rowData(sce)$CHR %in% c("MT", "X", "Y")), ]
# sce <- sce[which(!rowData(sce)$CHR %in% c("MT")), ]
set.seed(1011220)
clusters <- quickCluster(
  sce,
  min.size = 60,
  use.ranks = FALSE,
  BSPARAM = BiocSingular::IrlbaParam())
table(clusters)
sce <- computeSumFactors(
  sce,
  clusters = clusters,
  min.mean = 0.1)
summary(sizeFactors(sce))

xlim <- c(1, max(sce$sum) / 1e3)
ylim <- range(sizeFactors(sce))
par(mfrow = c(3, 3))
lapply(sort(unique(sce$plate_number)), function(p) {
  sce <- sce[, sce$plate_number == p]
  plot(
    x = sce$sum / 1e3,
    y = sizeFactors(sce),
    log = "xy",
    xlab = "Library size (thousands)",
    ylab = "Size factor",
    main = p,
    xlim = xlim,
    ylim = ylim,
    pch = 16)
})

sce <- computeSpikeFactors(sce, spikes = "ERCC", general.use = FALSE)
sce <- logNormCounts(sce)

var.out <- modelGeneVarWithSpikes(sce, "ERCC")
par(mfrow = c(1, 1))
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
fit <- metadata(var.out)
points(fit$mean, fit$var, col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes])
set.seed(1000)
sce <- denoisePCA(sce, technical=fit$trend, BSPARAM=BiocSingular::IrlbaParam())
ncol(reducedDim(sce, "PCA"))

set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50)
plotTSNE(sce, colour_by = "group")
plotTSNE(sce, colour_by = "cell_type_descriptor")
plotTSNE(sce, colour_by = "plate_number")
plotTSNE(sce, colour_by = "subsets_Mt_percent")

snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters <- cluster.out$membership
table(my.clusters)

sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster")

library(SingleR)
mouse_se <- MouseRNAseqData()
as.data.frame(colData(mouse_se)) %>%
  tabyl(label.main) %>%
  adorn_pct_formatting(1) %>%
  knitr::kable(
    caption = "Summary of main cell types in the `MouseRNAseqData` reference set.")
rmarkdown::paged_table(
  as.data.frame(colData(mouse_se)) %>%
    dplyr::count(label.main, label.fine) %>%
    dplyr::arrange(label.main))

pred_mouse_cell_main <- SingleR(
  test = sce,
  ref = mouse_se,
  labels = mouse_se$label.main,
  BPPARAM = SerialParam())
tabyl(data.frame(label.main = pred_mouse_cell_main$labels), label.main) %>%
  adorn_pct_formatting(digits = 1) %>%
  dplyr::arrange(desc(n)) %>%
  knitr::kable(
    caption = "Cell label assignments using the main labels of the `MouseRNAseqData` reference data.")
stopifnot(identical(rownames(pred_mouse_cell_main), colnames(sce)))
plotScoreHeatmap(
  pred_mouse_cell_main,
  annotation_col = data.frame(
    cluster = sce$cluster,
    sample = sce$cell_type_descriptor,
    row.names = rownames(pred_mouse_cell_main)))

pred_mouse_cell_fine <- SingleR(
  test = sce,
  ref = mouse_se,
  labels = mouse_se$label.fine,
  BPPARAM = SerialParam())
tabyl(data.frame(label.fine = pred_mouse_cell_fine$labels), label.fine) %>%
  adorn_pct_formatting(digits = 1) %>%
  dplyr::arrange(desc(n)) %>%
  knitr::kable(
    caption = "Cell label assignments using the fine labels of the `MouseRNAseqData` reference data.")
stopifnot(identical(rownames(pred_mouse_cell_main), colnames(sce)))
plotScoreHeatmap(
  pred_mouse_cell_fine,
  annotation_col = data.frame(
    cluster = sce$cluster,
    sample = sce$cell_type_descriptor,
    row.names = rownames(pred_mouse_cell_fine)))

genes_of_interest <- c(
  "Col4a1", "Col4a2", "Epas1", "Cdh5", "Ptprb", "Pecam1", "Vwf", "Itgb1",
  "Calcrl", "Plvap", "Tie1", "Cldn5", "Acvrl1", "Eng", "Kdr", "Kat7")
plotExpression(sce, genes_of_interest, x = "cell_type_descriptor", colour_by = "group")
data.frame(
  gene = genes_of_interest,
  median_expression = signif(rowMedians(as.matrix(logcounts(sce)[genes_of_interest, ])), 3)) %>%
  dplyr::arrange(desc(median_expression)) %>%
  knitr::kable()
p <- lapply(genes_of_interest, function(g) {
  plotTSNE(sce, colour_by = g) +
    ggtitle(g)
})
multiplot(plotlist = p, cols = 4)
