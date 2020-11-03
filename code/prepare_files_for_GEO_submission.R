# Prepare C075_Grant_Coultas data for GEO submission
# Peter Hickey
# 2020-11-03

library(here)

outdir <- here("tmp/GEO")
dir.create(outdir, recursive = TRUE)

# FASTQs -----------------------------------------------------------------------

dir.create(file.path(outdir, "FASTQ"))
# NOTE: These plate-level FASTQ files are created by code/scPipe.R
file.copy(
  from = here("extdata/NN158/scPipe/LC279/LC279.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN158/scPipe/LC279/LC279.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
for (plate in c("LC280", "LC294", "LC358", "LC392", "LC396", "LC398")) {
  message(plate)
  file.copy(
    from = here("extdata/NN176/scPipe", plate, paste0(plate, ".R1.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
  file.copy(
    from = here("extdata/NN176/scPipe", plate, paste0(plate, ".R2.fastq.gz")),
    to = file.path(outdir, "FASTQ"),
    recursive = FALSE,
    overwrite = FALSE)
}

# SCEs -------------------------------------------------------------------------

sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.scPipe.SCE.rds"))

# Update colData similar to that done in analysis/preprocessing.Rmd
sce$sample_type <- sce$sample_name
sce$sample_name <- NULL
sce$genotype <- vapply(
  strsplit(sce$cell_type_descriptor, " "),
  "[[",
  character(1),
  1)
sce$mouse <- gsub(
  " ",
  "",
  sub("(Control|Mutant) sample", "", sce$cell_type_descriptor))
mouse_sex <- read.csv(
  here("data", "sample_sheets", "mouse_sex.csv"),
  header = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE)
colnames(mouse_sex) <- c("mouse", "sex")
sce$sex <- dplyr::inner_join(
  dplyr::tibble(mouse = sce$mouse), mouse_sex) %>%
  dplyr::pull(sex)

# Gene counts
write.csv(
  x = as.data.frame(as.matrix(counts(sce))),
  file = gzfile(file.path(outdir, "SCE", "gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
write.csv(
  x = as.data.frame(as.matrix(counts(altExp(sce)))),
  file = gzfile(file.path(outdir, "SCE", "ERCC_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "sample_sheet.csv.gz")),
  row.names = TRUE)
