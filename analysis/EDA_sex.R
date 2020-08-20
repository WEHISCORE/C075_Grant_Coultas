library(scater)
library(scran)
library(here)
library(magrittr)

sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.scPipe.SCE.rds"))

mouse_sex <- read.csv(
  here("data", "sample_sheets", "mouse_sex.csv"),
  header = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE)
colnames(mouse_sex) <- c("mouse", "sex")

# Construct new_colData
new_colData <- DataFrame(
  plate_number = factor(sce$plate_number),
  well_position = sce$well_position,
  sample_type = factor(sce$sample_name),
  genotype = factor(
    vapply(strsplit(sce$cell_type_descriptor, " "), "[[", character(1), 1),
    c("Control", "Mutant")),
  mouse = factor(
    gsub(
      " ",
      "",
      sub("(Control|Mutant) sample", "", sce$cell_type_descriptor))),
  sequencing_run = factor(sce$sequencing_run),
  row.names = colnames(sce))
new_colData[["sex"]] <- dplyr::inner_join(
  dplyr::tibble(mouse = new_colData$mouse), mouse_sex) %>%
  dplyr::pull(sex) %>%
  factor()
# NOTE: The FACS data column names all start with 'B', 'V', 'Y', 'R', 'FSC' or
#       'SSC'.
facs_markers <- grep(
  "^B|^V|^R|^Y|^FSC|^SSC",
  colnames(colData(sce)),
  value = TRUE)
facs_data <- colData(sce)[, facs_markers]
# Assign new_colData to sce.
colData(sce) <- cbind(new_colData, facs_data)

# Extract rownames (Ensembl IDs) to use as key in database lookups.
ensembl <- rownames(sce)

# Pull out useful gene-based annotations from the Ensembl-based database.
library(EnsDb.Mmusculus.v79)
library(ensembldb)
# NOTE: These columns were customised for this project.
ensdb_columns <- c(
  "GENEBIOTYPE", "GENENAME", "GENESEQSTART", "GENESEQEND", "SEQNAME", "SYMBOL")
names(ensdb_columns) <- paste0("ENSEMBL.", ensdb_columns)
stopifnot(all(ensdb_columns %in% columns(EnsDb.Mmusculus.v79)))
ensdb_df <- DataFrame(
  lapply(ensdb_columns, function(column) {
    mapIds(
      x = EnsDb.Mmusculus.v79,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "GENEID",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)
# NOTE: Can't look up GENEID column with GENEID key, so have to add manually.
ensdb_df$ENSEMBL.GENEID <- ensembl

# NOTE: Mus.musculus combines org.Mm.eg.db and
#       TxDb.Mmusculus.UCSC.mm10.knownGene (as well as others) and therefore
#       uses entrez gene and RefSeq based data.
library(Mus.musculus)
# NOTE: These columns were customised for this project.
ncbi_columns <- c(
  # From TxDB: None required
  # From OrgDB
  "ALIAS", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL")
names(ncbi_columns) <- paste0("NCBI.", ncbi_columns)
stopifnot(all(ncbi_columns %in% columns(Mus.musculus)))
ncbi_df <- DataFrame(
  lapply(ncbi_columns, function(column) {
    mapIds(
      x = Mus.musculus,
      # NOTE: Need to remove gene version number prior to lookup.
      keys = gsub("\\.[0-9]+$", "", ensembl),
      keytype = "ENSEMBL",
      column = column,
      multiVals = "CharacterList")
  }),
  row.names = ensembl)

rowData(sce) <- cbind(ensdb_df, ncbi_df)

rownames(sce) <- uniquifyFeatureNames(
  ID = rownames(sce),
  # NOTE: An Ensembl ID may map to 0, 1, 2, 3, ... gene symbols.
  #       When there are multiple matches only the 1st match is used.
  names = sapply(rowData(sce)$ENSEMBL.SYMBOL, function(x) {
    if (length(x)) {
      x[[1]]
    } else {
      NA_character_
    }
  }))

ggcells(sce, aes(x = plate_number, y = log2(Xist + 1), col = sex), features = "Xist", exprs_values = "counts") +
  geom_violin(scale = "width", width = 1, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(~ sex) +
  cowplot::theme_minimal_hgrid() +
  guides(colour = FALSE)

ggcells(sce, aes(x = plate_number, y = log2(Tsix + 1), col = sex), features = "Tsix", exprs_values = "counts") +
  geom_violin(scale = "width", width = 1, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(~ sex) +
  cowplot::theme_minimal_hgrid() +
  guides(colour = FALSE)

ggcells(sce, aes(x = plate_number, y = log2(Ddx3y + 1), col = sex), features = "Ddx3y", exprs_values = "counts") +
  geom_violin(scale = "width", width = 1, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(~ sex) +
  cowplot::theme_minimal_hgrid() +
  guides(colour = FALSE)

ggcells(sce, aes(x = plate_number, y = log2(Eif2s3y + 1), col = sex), features = "Eif2s3y", exprs_values = "counts") +
  geom_violin(scale = "width", width = 1, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(~ sex) +
  cowplot::theme_minimal_hgrid() +
  guides(colour = FALSE)



ggcells(sce, aes(x = plate_number, y = log2(Ddx3x + 1), col = sex), features = "Ddx3x", exprs_values = "counts") +
  geom_violin(scale = "width", width = 1, alpha = 0.2) +
  ggbeeswarm::geom_quasirandom() +
  facet_grid(~ sex) +
  cowplot::theme_minimal_hgrid() +
  guides(colour = FALSE)

y <- rownames(rowData(sce)[any(rowData(sce)$ENSEMBL.SEQNAME == "Y"), ])
x <- rownames(rowData(sce)[any(rowData(sce)$ENSEMBL.SEQNAME == "X"), ])

sce <- logNormCounts(sce)

plotHeatmap(
  sce,
  features = c(
    x[rowMeans(logcounts(sce)[x, ]) > 1],
    y[rowMeans(logcounts(sce)[y, ]) > 0.05]),
  show_colnames = FALSE,
  center = TRUE,
  symmetric = TRUE,
  colour_columns_by = c("sex", "mouse"),
  order_columns_by = c("sex", "mouse"),
  cluster_rows = FALSE)

z <- findMarkers(sce, sce$sex, test.type = "binom", direction = "up")

plotHeatmap(
  sce,
  color = hcl.colors(101, "Blue-Red 3"),
  features = c(rownames(z$female)[1:10], rownames(z$male)[1:10]),
  show_colnames = FALSE,
  center = TRUE,
  symmetric = TRUE,
  colour_columns_by = c("sex","plate_number", "mouse"),
  order_columns_by = c("sex", "plate_number", "mouse"),
  cluster_rows = FALSE)
