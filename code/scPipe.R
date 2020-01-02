# Process NN158 (C075) with scPipe
# Peter Hickey
# 2020-01-02

# Setup ------------------------------------------------------------------------

library(scPipe)
library(Rsubread)
library(here)
library(readxl)
library(dplyr)
library(janitor)
library(stringr)
library(stringi)
library(tidyr)

options("mc.cores" = 8L)

# Construct NN158 sample sheet -------------------------------------------------

file_nn158 <- here(
  "data",
  "sample_sheets",
  "C078_and_C075_NN158_Wang_Cao_Zoe_Grant_SeqPrimer_layout_Aug19.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn158,
  sheet = "Sample & Index",
  skip = 2,
  n_max = 1)

# NOTE: FACS data in columns >= "K"
facs_data_idx <- seq(which(LETTERS == "K"), ncol(header_row))
header_row <- c(
  paste0(colnames(header_row[, -facs_data_idx]), header_row[1, -facs_data_idx]),
  unlist(header_row[1, facs_data_idx], use.names = FALSE))
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn158 <- read_excel(
  path = file_nn158,
  sheet = "Sample & Index",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)
# Restrict to Zoe's plates
sample_sheet_nn158 <- filter(
  sample_sheet_nn158,
  `Plate#` %in% c("LC279"))
# Tidy up names and empty rows/columns.
sample_sheet_nn158 <- bind_cols(
  clean_names(sample_sheet_nn158[, -facs_data_idx]),
  clean_names(sample_sheet_nn158[, facs_data_idx], case = "parsed"))
colnames(sample_sheet_nn158) <- gsub(
  "\\.\\.\\.",
  "",
  colnames(sample_sheet_nn158))
sample_sheet_nn158 <- remove_empty(sample_sheet_nn158)
colnames(sample_sheet_nn158) <- gsub(
  "\\_[0-9]+$",
  "",
  colnames(sample_sheet_nn158))

# Update FACS data (FACS data in columns >= "K")
facs_data_idx <- seq(which(LETTERS == "K"), ncol(sample_sheet_nn158))

# Filter out those without a cell index sequence, with no cell, or that were
# otherwise removed.
sample_sheet_nn158 <- sample_sheet_nn158 %>%
  filter(
    !is.na(rd1_index_cell_index_index_sequence_as_in_c_rt1_primer),
    !cell_type_descriptor %in% c("No Cell", "removed", "no cell"))

# Ensure FACS columns are stored as numeric (readxl sometimes fails, presumably
# to weird pattern of empty cells).
sample_sheet_nn158 <- sample_sheet_nn158 %>%
  mutate_at(facs_data_idx, as.numeric)

# Some final tidying.
sample_sheet_nn158 <- sample_sheet_nn158 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    sequencing_run = "NN158") %>%
  arrange(plate_number, well_position)

# Construct NN176 sample sheet -------------------------------------------------

file_nn176 <- here(
  "data",
  "sample_sheets",
  "C075_NN176_Zoe Grant_SeqPrimer layout_Dec19.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn176,
  sheet = "Sample & Index",
  skip = 2,
  n_max = 1)

# NOTE: FACS data in columns >= "J"
facs_data_idx <- seq(which(LETTERS == "J"), ncol(header_row))
header_row <- c(
  paste0(colnames(header_row[, -facs_data_idx]), header_row[1, -facs_data_idx]),
  unlist(header_row[1, facs_data_idx], use.names = FALSE))
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn176 <- read_excel(
  path = file_nn176,
  sheet = "Sample & Index",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)

# Tidy up names and empty rows/columns.
sample_sheet_nn176 <- bind_cols(
  clean_names(sample_sheet_nn176[, -facs_data_idx]),
  clean_names(sample_sheet_nn176[, facs_data_idx], case = "parsed"))
sample_sheet_nn176 <- remove_empty(sample_sheet_nn176)

# Filter out those without a cell index sequence, with no cell, or that were
# otherwise removed.
sample_sheet_nn176 <- sample_sheet_nn176 %>%
  filter(
    sample_name != "HRP control",
    !is.na(rd1_index_cell_index_index_sequence_as_in_c_rt1_primer) &
      !rd1_index_cell_index_index_sequence_as_in_c_rt1_primer %in%
      c("no primer", "No primer"),
    !cell_type_descriptor %in% c("No Cell", "empty", "no cell"))

# Ensure FACS columns are stored as numeric (readxl sometimes fails, presumably
# to weird pattern of empty cells).
sample_sheet_nn176 <- sample_sheet_nn176 %>%
  mutate_at(facs_data_idx, as.numeric)

# Some final tidying.
sample_sheet_nn176 <- sample_sheet_nn176 %>%
  mutate(
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    sequencing_run = "NN176") %>%
  arrange(plate_number, well_position)

# Construct final sample sheet -------------------------------------------------

sample_sheet <- bind_rows(sample_sheet_nn158, sample_sheet_nn176) %>%
  mutate(rowname = paste0(plate_number, "_", well_position)) %>%
  tibble::column_to_rownames("rowname") %>%
  DataFrame(., check.names = FALSE)

# NOTE: Check that there aren't any malformed well_positions (e.g., 'I19=A1').
stopifnot(!anyNA(sample_sheet$well_position))

# Key variables ----------------------------------------------------------------

plates <- unique(sample_sheet$plate_number)
names(plates) <- plates
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$plate_number,
  unique)
outdir <- here("data", "SCEs")
dir.create(outdir, recursive = TRUE)
extdir <- here("extdata", sequencing_runs, "scPipe", plates)
names(extdir) <- plates
sapply(extdir, dir.create, recursive = TRUE)
# NOTE: Only using first 7 nt of barcode.
read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
# NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
#       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
organism <- "mmusculus_gene_ensembl"
# NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
#       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
gene_id_type <- "ensembl_gene_id"

# Input files ------------------------------------------------------------------

# FASTQ files
r1_fq <- c(
  grep(
    pattern = "Undetermined",
    x = list.files(
      path = here("extdata", "NN158", "merged"),
      full.names = TRUE,
      pattern = glob2rx("*R1.fastq.gz")),
    invert = TRUE,
    value = TRUE),
  grep(
    pattern = "Undetermined",
    x = list.files(
      path = here("extdata", "NN176", "merged"),
      full.names = TRUE,
      pattern = glob2rx("*R1.fastq.gz")),
    invert = TRUE,
    value = TRUE))

r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))
tx_fq <- file.path(extdir, paste0(plates, ".R2.fastq.gz"))
names(tx_fq) <- plates
barcode_fq <- gsub("R2", "R1", tx_fq)

# NOTE: Concatenate FASTQs at the plate-level (this is safe because each plate
#       has a different pair of RPIs, one for the single-cells and one for the
#       20 cells).
mclapply(plates, function(plate) {
  message(plate)
  rpis <- sub(
    " ",
    "-",
    unique(
      sample_sheet$illumina_index_index_number_separate_index_read[sample_sheet$plate_number == plate]))
  cmd <- paste0(
    "cat ",
    paste0(
      sapply(
        rpis,
        function(x) grep(paste0(x, "_"), r1_fq, value = TRUE)),
      collapse = " "),
    " > ",
    barcode_fq[[plate]],
    "\n",
    "cat ",
    paste0(
      sapply(
        rpis,
        function(x) grep(paste0(x, "_"), r2_fq, value = TRUE)),
      collapse = " "),
    " > ",
    tx_fq[[plate]])
  system(cmd)
})

# Genome index
genome_index <- here("extdata", "GRCm38.p6", "GRCm38_with_ERCC")

# Genome annotation(s)
annofn <- c(
  here("extdata", "GRCm38.p6", "gencode.vM18.primary_assembly.annotation.gff3"),
  system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))

# Cell barcodes
bc_anno <- file.path(extdir, paste0(plates, ".barcode_annotation.csv"))
names(bc_anno) <- plates

for (plate in plates) {
  message(plate)
  tmp <- sample_sheet[sample_sheet$plate_number == plate, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    # NOTE: For some reason the primer name and sequence columns have been
    #       reversed in this sample sheet.
    # NOTE: Only using first 7 nt of barcode.
    barcode = strtrim(tmp$c_rt1_primer_name, 7),
    stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[plate]],
    quote = FALSE,
    row.names = FALSE)
}

# Output files -----------------------------------------------------------------

combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
names(combined_fq) <- names(tx_fq)
subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
exon_bam <- gsub("subread", "exon", subread_bam)

# FASTQ reformatting -----------------------------------------------------------

filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
mclapply(seq_along(tx_fq), function(i) {
  message(combined_fq[i])
  sc_trim_barcode(
    outfq = combined_fq[i],
    r1 = tx_fq[i],
    r2 = barcode_fq[i],
    read_structure = read_structure,
    filter_settings = filter_settings)
})

# Aligning reads to a reference genome -----------------------------------------

Rsubread::align(
  index = genome_index,
  readfile1 = combined_fq,
  output_file = subread_bam,
  nthreads = 24)

# Assigning reads to annotated exons -------------------------------------------

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
mclapply(seq_along(subread_bam), function(i) {
  message(i)
  sc_exon_mapping(
    inbam = subread_bam[i],
    outbam = exon_bam[i],
    annofn = annofn,
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
})

# De-multiplexing data ---------------------------------------------------------

max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
mclapply(seq_along(exon_bam), function(i) {
  message(i)
  sc_demultiplex(
    inbam = exon_bam[i],
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)
})

# Gene counting ----------------------------------------------------------------

UMI_cor <- 1
gene_fl <- FALSE
mclapply(seq_along(bc_anno), function(i) {
  message(i)
  sc_gene_counting(
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
})

# Create and save SingleCellExperiment -----------------------------------------

list_of_sce <- lapply(plates, function(plate) {
  create_sce_by_dir(
    datadir = extdir[[plate]],
    organism = organism,
    gene_id_type = gene_id_type,
    pheno_data = sample_sheet[sample_sheet$plate_number == plate, ],
    # NOTE: Create the report separately for more fine-grained control.
    report = FALSE)
})

source(here("analysis", "helper_functions.R"))
sce <- Reduce(function(x, y) .combine(x, y, rowData_by = NULL), list_of_sce)
assay(sce, withDimnames = FALSE) <- as(
  assay(sce, withDimnames = FALSE),
  "dgCMatrix")
sce <- splitAltExps(sce, ifelse(isSpike(sce), "ERCC", "Endogenous"))
sce <- clearSpikes(sce)
saveRDS(
  sce,
  file.path(outdir, "C075_Grant_Coultas.scPipe.SCE.rds"),
  compress = "xz")

# Create QC report -------------------------------------------------------------

# NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
dir.create(here("output", "scPipe"), recursive = TRUE)
# NOTE: Tends to crap itself if using mclapply().
lapply(plates, function(plate) {
  try(create_report(
    sample_name = plate,
    outdir = extdir[[plate]],
    r1 = tx_fq[[plate]],
    r2 = barcode_fq[[plate]],
    outfq = combined_fq[[plate]],
    read_structure = read_structure,
    filter_settings = filter_settings,
    align_bam = subread_bam[[plate]],
    genome_index = genome_index,
    map_bam = exon_bam[[plate]],
    exon_anno = annofn,
    stnd = stnd,
    fix_chr = fix_chr,
    barcode_anno = bc_anno[[plate]],
    max_mis = max_mis,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl,
    organism = organism,
    gene_id_type = gene_id_type))

  # NOTE: Workaround bug in create_report() and stop output after 'Data summary'
  #       section.
  tmp <- readLines(file.path(extdir[[plate]], "report.Rmd"))
  tmp <- c(tmp[1:161], "knitr::knit_exit()", tmp[162:length(tmp)])
  writeLines(tmp, file.path(extdir[[plate]], "report.Rmd"))
  knitr::wrap_rmd(
    file = file.path(extdir[[plate]], "report.Rmd"),
    width = 120,
    backup = NULL)
  rmarkdown::render(
    input = file.path(extdir[[plate]], "report.Rmd"),
    output_file = file.path(extdir[[plate]], "report.html"),
    knit_root_dir = ".")

  # NOTE: Copy the QC report to the repository.
  file.copy(
    from = file.path(extdir[[plate]], "report.nb.html"),
    to = here(
      "output",
      "scPipe",
      paste0(plate, ".scPipe_QC_report.nb.html")),
    overwrite = TRUE)
})
