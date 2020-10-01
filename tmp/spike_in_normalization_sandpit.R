# HBO1 causes global expression changes and Anne would like to quantify these.
# Here, I experimenting with ideas in
# https://osca.bioconductor.org/normalization.html#spike-norm to see if this
# can help us achieve that.

library(SingleCellExperiment)
library(here)
library(scater)
library(scran)

sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.preprocessed.SCE.rds"))

# Some useful colours
plate_number_colours <- setNames(
  unique(sce$plate_number_colours),
  unique(names(sce$plate_number_colours)))
genotype_colours <- setNames(
  unique(sce$genotype_colours),
  unique(names(sce$genotype_colours)))
mouse_colours <- setNames(
  unique(sce$mouse_colours),
  unique(names(sce$mouse_colours)))
sex_colours <- setNames(
  unique(sce$sex_colours),
  unique(names(sce$sex_colours)))
cluster_colours <- setNames(
  unique(sce$cluster_colours),
  unique(names(sce$cluster_colours)))

# EDA --------------------------------------------------------------------------

plotColData(
  sce,
  x = "genotype",
  y = "sum",
  other_fields = c("plate_number"),
  colour_by = "genotype") +
  scale_y_log10(limits = c(1, NA)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm")) +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  stat_summary(fun.data = "mean_cl_boot")

plotColData(
  sce,
  x = "genotype",
  y = I(colSums(counts(altExp(sce, "ERCC")))),
  other_fields = c("plate_number"),
  colour_by = "genotype") +
  scale_y_log10(limits = c(1, NA)) +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm")) +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  stat_summary(fun.data = "mean_cl_boot")

plotColData(
  sce,
  x = "genotype",
  y = I(colSums(counts(altExp(sce, "ERCC"))) / sce$sum),
  other_fields = c("plate_number"),
  colour_by = "genotype") +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 0.05, lty = 2, col = "red")

# Normalization by spike-ins ---------------------------------------------------

sce <- computeSpikeFactors(sce, "ERCC")

plotColData(
  sce,
  x = "genotype",
  y = I(sizeFactors(sce)),
  other_fields = c("plate_number"),
  colour_by = "genotype") +
  scale_y_log10() +
  annotation_logticks(
    sides = "l",
    short = unit(0.03, "cm"),
    mid = unit(0.06, "cm"),
    long = unit(0.09, "cm")) +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  stat_summary(fun.data = "mean_cl_boot") +
  geom_hline(yintercept = 1, lty = 2, col = "red")

clusters <- quickCluster(sce)
clusters2 <- quickCluster(sce, block = sce$plate_number)
to.plot <- as.data.frame(
  cbind(
    data.frame(
      DeconvFactor=calculateSumFactors(sce, clusters = clusters2, min.mean = 0.1),
      SpikeFactor=sizeFactors(sce)),
    colData(sce)))

ggplot(to.plot, aes(x = DeconvFactor, y = SpikeFactor, colour = plate_number)) +
  geom_point() +
  facet_grid(~genotype) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_colour_manual(values = plate_number_colours)

ggplot(to.plot, aes(x = DeconvFactor, y = SpikeFactor, colour = plate_number)) +
  geom_point() +
  facet_grid(sex~genotype) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_colour_manual(values = plate_number_colours)

ggplot(to.plot,
  aes(x = DeconvFactor, y = SpikeFactor, colour = genotype)) +
  geom_point() +
  facet_grid(~plate_number) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_colour_manual(values = genotype_colours)

ggplot(
  as.data.frame(to.plot),
  aes(x = DeconvFactor, y = SpikeFactor, colour = genotype)) +
  geom_point() +
  facet_grid(genotype~plate_number) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  scale_colour_manual(values = genotype_colours)

sce.deconv <- sce
sizeFactors(sce.deconv) <- to.plot$DeconvFactor
sce.spike <- sce
sizeFactors(sce) <- to.plot$SpikeFactor
sce.deconv <- logNormCounts(sce.deconv)
sce.spike <- logNormCounts(sce.spike)

a <- list(
  plotExpression(sce.deconv, x="genotype",
                 other_fields = "plate_number", features="Arap2") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After deconvolution") +
    facet_grid(~plate_number),
  plotExpression(sce.spike, x="genotype",
                 other_fields="plate_number", features="Arap2") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After spike-in normalization") +
    facet_grid(~plate_number))

sce.deconv <- batchelor::multiBatchNorm(sce.deconv, batch = sce$plate_number)
sce.spike <- batchelor::multiBatchNorm(sce.spike, batch = sce$plate_number)

b <- list(
  plotExpression(sce.deconv, x="genotype",
                 other_fields = "plate_number", features="Arap2") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After deconvolution") +
    facet_grid(~plate_number),
  plotExpression(sce.spike, x="genotype",
                 other_fields="plate_number", features="Arap2") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After spike-in normalization") +
    facet_grid(~plate_number))

cowplot::plot_grid(plotlist= c(a, b), ncol = 2)

sce <- logNormCounts(sce)

# Merging ----------------------------------------------------------------------

sce <- readRDS("data/SCEs/C075_Grant_Coultas.preprocessed.SCE.rds")
clusters <- quickCluster(sce)
clusters2 <- quickCluster(sce, block = sce$plate_number)
sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1)

sce$batch <- sce$plate_number
# NOTE: Can reuse plate colours
batch_colours <- plate_number_colours

dec <- modelGeneVarWithSpikes(sce, "ERCC", block = sce$batch)

xlim <- c(0, max(sapply(dec$per.block, function(x) max(x$mean))))
ylim <- c(0, max(sapply(dec$per.block, function(x) max(x$total))))
par(mfrow = c(2, 3))
blocked_stats <- dec$per.block
for (i in colnames(blocked_stats)) {
  current <- blocked_stats[[i]]
  plot(
    current$mean,
    current$total,
    main = i,
    pch = 16,
    cex = 0.5,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression",
    xlim = xlim,
    ylim = ylim)
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col = "red", pch = 16)
  curve(curfit$trend(x), col = batch_colours[[i]], add = TRUE, lwd = 2)
}

hvg <- getTopHVGs(dec, n = 2000)

set.seed(2719)
sce <- denoisePCA(sce, technical = dec, subset.row = hvg)
dim(reducedDim(sce))
attr(reducedDim(sce), "percentVar")

library(batchelor)
mnn_out <- fastMNN(
  multiBatchNorm(sce, batch = factor(sce$plate_number)),
  batch = factor(sce$plate_number),
  cos.norm = FALSE,
  d = ncol(reducedDim(sce, "PCA")),
  auto.merge = TRUE,
  # merge.order = list(
  #   # list("LC358", "LC294", "LC392", "LC280"),
  #   list("LC392", "LC280"),
  #   list("LC396", "LC398"),
  #   list("LC279")),
  subset.row = hvg)
var_loss_auto <- metadata(mnn_out)$merge.info$lost.var
rownames(var_loss_auto) <- paste0("merge", 1:nrow(var_loss_auto))
knitr::kable(
  100 * var_loss_auto,
  digits = 1)

reducedDim(sce, "corrected") <- reducedDim(mnn_out)

set.seed(572968)
library(uwot)
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
  ggtitle("MNN-corrected") +
  xlab("Dimension 1") +
  ylab("Dimension 2")

# A different EDA --------------------------------------------------------------

set.seed(666)
sce_deconv <- computeSumFactors(sce, cluster = quickCluster(sce), min.mean = 0.1)
sce_spike <- computeSpikeFactors(sce, "ERCC")
cce_mbn_spike <- batchelor::multiBatchNorm(
  sce_spike,
  batch = sce_spike$plate_number,
  min.mean = 0.1)
sce_mbn_deconv <- batchelor::multiBatchNorm(
  sce_deconv,
  batch = sce_deconv$plate_number,
  min.mean = 0.1)
sce_spike_mbn <- sce_spike
sizeFactors(sce_spike_mbn) <- sizeFactors(
  batchelor::multiBatchNorm(
    altExp(sce_spike),
    batch = sce_spike$plate_number,
    min.mean = 0.1))

sce_spike_mynorm <- sce_spike
sizeFactors(sce_spike_mynorm) <- sizeFactors(sce_spike_mynorm) /
  tapply(sizeFactors(sce_spike_mynorm), sce_spike_mynorm$plate_number, median)[sce_spike_mynorm$plate_number]

sf_df <- data.frame(
  deconv = sizeFactors(sce_deconv),
  spike = sizeFactors(sce_spike),
  mbn_deconv = sizeFactors(sce_mbn_deconv),
  mbn_spike = sizeFactors(sce_mbn_spike),
  spike_mnb = sizeFactors(sce_spike_mbn),
  spike_mynorm = sizeFactors(sce_spike_mynorm))
my_line <- function(x, y, ...) {
  points(x, y, ...)
  abline(a = 0, b = 1, col = "red")
}
pairs(
  sf_df,
  col = sce$plate_number_colours,
  log = "xy",
  xlim = range(sf_df),
  ylim = range(sf_df),
  lower.panel = my_line,
  upper.panel = my_line)

to.plot <- data.frame(
  DeconvFactor=sizeFactors(sce_deconv),
  SpikeFactor=sizeFactors(sce_spike),
  MBNDeconvFactor=sizeFactors(sce_mbn_deconv),
  MBNSpikeFactor=sizeFactors(sce_mbn_spike),
  SpikeFactorMBN=sizeFactors(sce_spike_mbn),
  SpikeFactorMyNorm=sizeFactors(sce_spike_mynorm),
  Genotype=sce$genotype,
  Plate=sce$plate_number)

q <- ggplot(to.plot, aes(x=DeconvFactor, y=SpikeFactor, color=Genotype)) +
  geom_point() + facet_wrap(~Plate) + scale_x_log10() +
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red") +
  scale_colour_manual(values = genotype_colours)

qq <- ggplot(to.plot, aes(x = Genotype, y=SpikeFactor, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()

qqq <- ggplot(to.plot, aes(x = Genotype, y=DeconvFactor, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()

p <- plotColData(
  sce,
  x = "genotype",
  y = "altexps_ERCC_percent",
  other_fields = "plate_number",
  colour_by = "genotype") +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  # ylim(0, max(sce$altexps_ERCC_percent)) +
  ggtitle("C075") +
  ylab("ERCC percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pp <- plotColData(
  sce,
  x = "genotype",
  y = "sum",
  other_fields = "plate_number",
  colour_by = "genotype") +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  scale_y_log10() +
  ggtitle("C075") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ppp <- plotColData(
  sce,
  x = "genotype",
  y = "altexps_ERCC_sum",
  other_fields = "plate_number",
  colour_by = "genotype") +
  scale_fill_manual(values = genotype_colours) +
  facet_grid(~plate_number) +
  guides(fill = FALSE) +
  scale_y_log10() +
  ggtitle("C075") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

qqqq <- ggplot(to.plot, aes(x = Genotype, y=MBNSpikeFactor, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()
qqqqq <- ggplot(to.plot, aes(x = Genotype, y=MBNDeconvFactor, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()

qqqqqq <- ggplot(to.plot, aes(x = Genotype, y=SpikeFactorMBN, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()

qqqqqqq <- ggplot(to.plot, aes(x = Genotype, y=SpikeFactorMyNorm, colour=Genotype)) +
  geom_violin() + facet_grid(~Plate) + scale_colour_manual(values = genotype_colours) + guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10()

p + qq + qqq +
  pp + qqqq + qqqqq +
  ppp + qqqqqq + qqqqqqq +
  plot_layout(ncol = 3)
