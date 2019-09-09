library(here)
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.scPipe.SCE.rds"))

library(scater)
library(dplyr)
library(janitor)
sce <- calculateQCMetrics(sce)

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
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = which(is_mito)))

sce$sample <- paste0(
  sce$sample_name,
  ".",
  ifelse(grepl("Control", sce$cell_type_descriptor), "control", "mutant"))
multiplot(
  plotColData(sce, "log10_total_counts", x = "sample", colour_by = "cell_type_descriptor") +
    ylim(0, NA) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  plotColData(sce, "log10_total_features_by_counts", x = "sample", colour_by = "cell_type_descriptor") +
    ylim(0, NA) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  plotColData(sce, "pct_counts_ERCC", x = "sample", colour_by = "cell_type_descriptor") +
    ylim(0, 100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  plotColData(sce, "pct_counts_Mt", x = "sample", colour_by = "cell_type_descriptor") +
    ylim(0, 100) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  cols = 2)

pdf(here("tmp", "quickQC.pdf"), width = 12, height = 6)
plotPlatePosition(sce, sce$well_position, shape_by = "sample", point_size = 6, colour_by = "log10_total_counts") +
  scale_shape_manual(values = c("20 cells.control" = 15, "20 cells.mutant" = 16, "single cell.control" = 17, "single cell.mutant" = 18))
plotPlatePosition(sce, sce$well_position, shape_by = "sample", point_size = 6, colour_by = "log10_total_features_by_counts") +
  scale_shape_manual(values = c("20 cells.control" = 15, "20 cells.mutant" = 16, "single cell.control" = 17, "single cell.mutant" = 18))
plotPlatePosition(sce, sce$well_position, shape_by = "sample", point_size = 6, colour_by = "pct_counts_ERCC") +
  scale_shape_manual(values = c("20 cells.control" = 15, "20 cells.mutant" = 16, "single cell.control" = 17, "single cell.mutant" = 18))
plotPlatePosition(sce, sce$well_position, shape_by = "sample", point_size = 6, colour_by = "pct_counts_Mt") +
  scale_shape_manual(values = c("20 cells.control" = 15, "20 cells.mutant" = 16, "single cell.control" = 17, "single cell.mutant" = 18))
dev.off()

# Check expression of mutant gene (Kat7)
plotExpression(sce, "Kat7", x = "sample", exprs_values = "counts")
table(sce$sample, counts(sce)["Kat7", ] > 0)

sce <- sce[, sce$sample_name == "single cell"]
sce$sample <- factor(sce$sample)
sce$sample_name <- factor(sce$sample_name)

libsize_drop <- isOutlier(
  metric = sce$total_counts,
  nmads = 3,
  type = "lower",
  log = TRUE)
feature_drop <- isOutlier(
  metric = sce$total_features_by_counts,
  nmads = 3,
  type = "lower",
  log = TRUE)
spike_drop <- isOutlier(
  metric = sce$pct_counts_ERCC,
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
  "pre-QC" = table(sce_pre_QC_outlier_removal$sample),
  "post-QC" = table(sce$sample))

plotHighestExprs(sce, n = 50)
plotHighestExprs(sce[!grepl("^ERCC", rownames(sce))], n = 50)

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
  # NOTE: Data are small enough to use exact PCA.
  BSPARAM = BiocSingular::ExactParam())
table(clusters)
sce <- computeSumFactors(
  sce,
  clusters = clusters,
  min.mean = 0.1)
summary(sizeFactors(sce))

plot(
  x = sce$total_counts / 1e3,
  y = sizeFactors(sce),
  log = "xy",
  xlab = "Library size (thousands)",
  ylab = "Size factor",
  pch = 16)

sce <- computeSpikeFactors(sce, type = "ERCC", general.use = FALSE)
sce <- normalize(sce)

var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.4))
var.out <- decomposeVar(sce, var.fit)
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
points(var.out$mean[isSpike(sce)], var.out$total[isSpike(sce)], col="red", pch=16)
curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out)[chosen.genes])
set.seed(1000)
sce <- denoisePCA(sce, technical=var.fit$trend, BSPARAM=BiocSingular::ExactParam())
ncol(reducedDim(sce, "PCA"))

set.seed(1000)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=50)
plotTSNE(sce, colour_by = "sample")

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
    sample = sce$sample,
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
    sample = sce$sample,
    row.names = rownames(pred_mouse_cell_fine)))

genes_of_interest <- c(
  "Col4a1", "Col4a2", "Epas1", "Cdh5", "Ptprb", "Pecam1", "Vwf", "Itgb1",
  "Calcrl", "Plvap", "Tie1", "Cldn5", "Acvrl1", "Eng", "Kdr", "Kat7")
plotExpression(sce, genes_of_interest, x = "sample")
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
