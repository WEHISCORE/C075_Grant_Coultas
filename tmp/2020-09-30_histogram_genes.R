library(scater)
library(here)
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.annotated.SCE.rds"))

# Violin plots -----------------------------------------------------------------

label_colours <- setNames(
  unique(sce$label_colours),
  unique(names(sce$label_colours)))

histogram_genes <- readxl::read_excel(
  here("data/zoes_marker_genes/histogram_genes.xlsx"))

pdf(here("tmp/2020-09-30_violin_plots.pdf"), width = 4, height = 8)
for (g in histogram_genes$Gene) {
  message(g)
  p <- plotExpression(
    sce,
    g,
    x =  "label",
    other_fields = c("genotype", "plate_number"),
    colour_by = "label",
    show_median = TRUE) +
    scale_fill_manual(values = label_colours) +
    facet_grid(plate_number ~ genotype) +
    ggtitle(g) +
    cowplot::theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = FALSE)
  print(p)
}
dev.off()

# Heatmaps ---------------------------------------------------------------------

# Define colData columns to be kept for pseudo-bulked data.
columns_to_keep <- c("mouse", "plate_number", "genotype", "sex")

summed <- aggregateAcrossCells(
  sce,
  id = colData(sce)[, c("label", "sample", columns_to_keep)],
  coldata_merge = FALSE)

summed <- logNormCounts(summed)
colnames(summed) <- paste0(summed$label, ".", summed$sample)

o <- order(summed$label, summed$genotype)
features <- histogram_genes$Gene

mat <- logcounts(summed)[features, ]
mat <- mat - rowMeans(mat)
zlim <- c(-3, 3)
mat[mat < zlim[1]] <- zlim[1]
mat[mat > zlim[2]] <- zlim[2]
p1 <- pheatmap(
  mat[, o],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  annotation_col = data.frame(
    label = summed$label[o],
    genotype = summed$genotype[o],
    row.names = colnames(summed)[o]),
  annotation_row = data.frame(
    type = histogram_genes$`Endothelial cell type`,
    row.names = histogram_genes$Gene),
  annotation_colors = list(
    label = label_colours,
    genotype = genotype_colours,
    plate_number = plate_number_colours,
    sex = sex_colours),
  show_colnames = FALSE,
  main = "row-normalized",
  silent = TRUE)

mat <- logcounts(summed)[features, ]
for (label in levels(sce$label)) {
  j <- grep(label, colnames(mat))
  mat[, j] <- mat[, j] - rowMeans(mat[, j])
}
zlim <- c(-2, 2)
mat[mat < zlim[1]] <- zlim[1]
mat[mat > zlim[2]] <- zlim[2]
o <- order(summed$label, summed$genotype, summed$sample)
p2 <- pheatmap(
  mat[, o],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  annotation_col = data.frame(
    label = summed$label[o],
    genotype = summed$genotype[o],
    row.names = colnames(summed)[o]),
  annotation_row = data.frame(
    type = histogram_genes$`Endothelial cell type`,
    row.names = histogram_genes$Gene),
  annotation_colors = list(
    label = label_colours,
    genotype = genotype_colours,
    plate_number = plate_number_colours,
    sex = sex_colours),
  show_colnames = FALSE,
  main = "row-normalized by label",
  silent = TRUE)

gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2)

# Heatmap of LFCs (as suggested by Luyi)
list_of_lfcs <- lapply(de_results, function(res) {
  res[features, "logFC"]
})
mat <- do.call(cbind, list_of_lfcs[levels(summed$label)])
rownames(mat) <- features
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  annotation_row = data.frame(
    type = histogram_genes$`Endothelial cell type`,
    row.names = histogram_genes$Gene),
  annotation_colors = list(
    label = label_colours,
    genotype = genotype_colours,
    plate_number = plate_number_colours,
    sex = sex_colours))

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hcl.colors(101, "Blue-Red 3"),
  annotation_row = data.frame(
    type = histogram_genes$`Endothelial cell type`,
    row.names = histogram_genes$Gene),
  annotation_colors = list(
    label = label_colours,
    genotype = genotype_colours,
    plate_number = plate_number_colours,
    sex = sex_colours))
