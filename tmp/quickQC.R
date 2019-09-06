library(here)
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.scPipe.SCE.rds"))

library(scater)
library(dplyr)
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
