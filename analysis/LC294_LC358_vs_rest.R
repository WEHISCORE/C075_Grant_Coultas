# After preprocessing ----------------------------------------------------------

library(scater)
library(scran)
library(here)

source(here("code/helper_functions.R"))
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.preprocessed.SCE.rds"))

plate_number_colours <- setNames(
  RColorBrewer::brewer.pal(nlevels(sce$plate_number), "Set2"),
  levels(sce$plate_number))

sce$groups <- ifelse(
  sce$plate_number %in% c("LC294", "LC358"),
  "LC294_LC358",
  "not_LC294_LC358")

markers <- findMarkers(sce, groups = sce$groups)
write.csv(as.data.frame(
  markers$LC294_LC358),
  here("tmp", "LC294_LC358_vs_rest.DEGs.csv"))

summary(markers$LC294_LC358$FDR < 0.05)

pdf(here("tmp/LC294_LC358_vs_rest.DEGs.pdf"), height = 12, width = 8)
.plotHeatmap(
  sce,
  features = rownames(markers$LC294_LC358)[1:100],
  center = TRUE,
  colour_columns_by = "plate_number",
  annotation_colors = list(plate_number = plate_number_colours),
  zlim = c(-3, 3),
  symmetric = TRUE,
  order_columns_by = "plate_number")
dev.off()

markers_wilcox_up <- findMarkers(
  sce,
  groups = sce$groups,
  test = "wilcox",
  direction = "up")
write.csv(as.data.frame(
  markers_wilcox_up$LC294_LC358),
  here("tmp", "LC294_LC358_vs_rest.DEGs.wilcox_up.csv"))
write.csv(as.data.frame(
  markers_wilcox_up$not_LC294_LC358),
  here("tmp", "rest_vs_LC294_LC358.DEGs.wilcox_up.csv"))

pdf(here("tmp/LC294_LC358_vs_rest.DEGs.wilcox_up.pdf"), height = 12, width = 8)
.plotHeatmap(
  sce,
  features = rownames(markers_wilcox_up$LC294_LC358)[1:100],
  center = TRUE,
  colour_columns_by = "plate_number",
  annotation_colors = list(plate_number = plate_number_colours),
  zlim = c(-3, 3),
  symmetric = TRUE,
  order_columns_by = "plate_number")
dev.off()

pdf(here("tmp/rest_vs_LC294_LC358.DEGs.wilcox_up.pdf"), height = 12, width = 8)
.plotHeatmap(
  sce,
  features = rownames(markers_wilcox_up$not_LC294_LC358)[1:100],
  center = TRUE,
  colour_columns_by = "plate_number",
  annotation_colors = list(plate_number = plate_number_colours),
  zlim = c(-3, 3),
  symmetric = TRUE,
  order_columns_by = "plate_number")
dev.off()

markers_binom_up <- findMarkers(
  sce,
  groups = sce$groups,
  test = "binom",
  direction = "up")
write.csv(as.data.frame(
  markers_binom_up$LC294_LC358),
  here("tmp", "LC294_LC358_vs_rest.DEGs.binom_up.csv"))

pdf(here("tmp/LC294_LC358_vs_rest.DEGs.binom_up.pdf"), height = 12, width = 8)
.plotHeatmap(
  sce,
  features = rownames(markers_binom_up$LC294_LC358)[1:100],
  center = FALSE,
  colour_columns_by = "plate_number",
  annotation_colors = list(plate_number = plate_number_colours),
  symmetric = FALSE,
  order_columns_by = "plate_number",
  zlim = c(0, 5))
dev.off()

# Removing cell cycle genes ----------------------------------------------------

# Adapted from biomaRt vignette
library(biomaRt)

ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
go <- c("GO:0007049")
ens <- getBM(
  attributes = "ensembl_gene_id",
  filters = "go",
  values = list(go),
  mart = ensembl)

sym <- rownames(sce)[rowData(sce)$ENSEMBL %in% ens$ensembl_gene_id]
sym <- union(sym, grep("^Hist", rownames(sce), value = TRUE))

hvg <- getTopHVGs(dec[setdiff(rownames(dec), sym), ], n = 1000)

plotExpression(object = sce, features = hvg[1:10])

sce <- denoisePCA(sce, technical = dec, subset.row = hvg)

library(batchelor)
mnn_out <- fastMNN(
  multiBatchNorm(sce, batch = sce$plate_number),
  batch = sce$plate_number,
  cos.norm = FALSE,
  d = ncol(reducedDim(sce, "PCA")),
  # auto.merge = TRUE,
  merge.order = list(
    list("LC358", "LC294", "LC392", "LC280"),
    list("LC396", "LC398"),
    list("LC279")),
  subset.row = hvg)
metadata(mnn_out)$merge.info
metadata(mnn_out)$merge.info$lost.var

reducedDim(sce, "corrected") <- reducedDim(mnn_out)

set.seed(572968)
uncorrected <- runUMAP(sce, dimred = "PCA")
uncorrected_umap <- cbind(
  data.frame(
    x = reducedDim(uncorrected, "UMAP")[, 1],
    y = reducedDim(uncorrected, "UMAP")[, 2]),
  as.data.frame(colData(uncorrected)))

ggplot(aes(x = x, y = y), data = uncorrected_umap) +
  geom_point(aes(colour = plate_number), size = 1) +
  scale_colour_manual(values = plate_number_colours) +
  cowplot::theme_cowplot(font_size = 10) +
  ggtitle("Original") +
  xlab("Dimension 1") +
  ylab("Dimension 2")

set.seed(55788)
sce <- runUMAP(sce, dimred = "corrected")
corrected_umap <- cbind(
  data.frame(
    x = reducedDim(sce, "UMAP")[, 1],
    y = reducedDim(sce, "UMAP")[, 2]),
  as.data.frame(colData(sce)))

ggplot(aes(x = x, y = y), data = corrected_umap) +
  geom_point(aes(colour = plate_number), size = 1) +
  scale_colour_manual(values = plate_number_colours) +
  cowplot::theme_cowplot(font_size = 10) +
  ggtitle("Original") +
  xlab("Dimension 1") +
  ylab("Dimension 2")

