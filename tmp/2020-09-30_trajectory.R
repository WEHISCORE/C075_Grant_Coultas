# Quick pseudotime/trajectory analysis adapted from https://osca.bioconductor.org/trajectory-analysis.html
# Peter Hickey
# 2020-09-30

# NOTE: This didn't prove very useful.

library(here)
sce <- readRDS(here("data/SCEs/C075_Grant_Coultas.annotated.SCE.rds"))

# TSCAN ------------------------------------------------------------------------

library(scater)
by.cluster <- aggregateAcrossCells(sce, ids = sce$cluster)
centroids <- reducedDim(by.cluster, "corrected")

dmat <- dist(centroids)
dmat <- as.matrix(dmat)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)

set.seed(1000)
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "UMAP")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)

plotUMAP(sce, colour_by="cluster", text_by = "cluster") +
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))

# TODO: for the love of god, we definitely need to move this into a function!
.map2edges <- function(points, center, edge.ends, previous) {
  all.distances <- list()
  all.pseudo <- list()
  edge.len <- list()

  # Computing distance of each point from each edge.
  # Edges defined from 'center' to 'edge.ends'.
  for (i in rownames(edge.ends)) {
    edge.end <- edge.ends[i,]
    delta <- center - edge.end
    max.d <- sqrt(sum(delta^2))
    delta <- delta/max.d

    centered <- t(t(points) - center)
    proj <- as.numeric(centered %*% delta)
    proj <- pmax(0, pmin(proj, max.d))
    mapped <- outer(proj, delta)

    dist <- sqrt(rowSums((centered - mapped)^2))
    all.distances[[i]] <- dist
    all.pseudo[[i]] <- proj
    edge.len[[i]] <- max.d
  }

  all.distances <- do.call(cbind, all.distances)
  all.pseudo <- do.call(cbind, all.pseudo)
  chosen <- colnames(all.distances)[max.col(-all.distances)]

  # Flipping the distance of points to the previous node,
  # in order to enforce a directional pseudo-time.
  dist.previous <- 0
  if (!is.na(previous)) {
    on.previous <- chosen==previous
    dist.previous <- edge.len[[previous]]
    previous.proj <- dist.previous - all.pseudo[on.previous,previous,drop=FALSE]

    if (all(on.previous)) {
      return(list(dist=dist.previous, pseudo=list(previous.proj)))
    }
  }

  # Filling out the branches, where points are NA for a branch's
  # pseudo-time if they were assigned to another branch.
  output <- list()
  for (leftover in setdiff(rownames(edge.ends), previous)) {
    empty <- rep(NA_real_, nrow(points))
    if (!is.na(previous)) {
      empty[on.previous] <- previous.proj
    }
    current <- chosen==leftover
    empty[current] <- all.pseudo[current,leftover]
    output[[leftover]] <- empty
  }

  list(dist=dist.previous, pseudo=output)
}

originals <- reducedDim(sce, "corrected")
cluster <- sce$cluster
starting.cluster <- names(igraph::V(mst)[igraph::degree(mst)==1])[1]
collated <- list()

latest <- starting.cluster
parents <- NA_character_
progress <- list(rep(NA_real_, length(cluster)))
cumulative <- 0

while (length(latest)) {
  new.latest <- new.parents <- character(0)
  new.progress <- list()
  new.cumulative <- numeric(0)

  for (i in seq_along(latest)) {
    curnode <- latest[i]
    all.neighbors <- names(igraph::adjacent_vertices(mst, curnode, mode="all")[[1]])
    in.cluster <- cluster==curnode

    mapped <- .map2edges(originals[in.cluster,,drop=FALSE], center=centroids[curnode,],
                         edge.ends=centroids[all.neighbors,,drop=FALSE], previous=parents[i])
    edge.len <- mapped$dist
    pseudo <- mapped$pseudo

    collected.progress <- list()
    for (j in seq_along(pseudo)) {
      sofar <- progress[[i]] # yes, using 'i' here.
      sofar[in.cluster] <- pseudo[[j]] + cumulative[i]
      collected.progress[[j]] <- sofar
    }

    all.children <- setdiff(all.neighbors, parents[i])
    if (length(all.children)==0) {
      collated[[curnode]] <- collected.progress[[1]]
    } else {
      new.latest <- c(new.latest, all.children)
      new.parents <- c(new.parents, rep(curnode, length(all.children)))
      new.progress <- c(new.progress, collected.progress)
      new.cumulative <- c(new.cumulative, rep(cumulative[i] + edge.len, length(all.children)))
    }
  }

  latest <- new.latest
  parents <- new.parents
  progress <- new.progress
  cumulative <- new.cumulative
}
tscan.pseudo <- do.call(cbind, collated)

plotUMAP(sce, colour_by=I(rowMeans(tscan.pseudo, na.rm=TRUE)), text_by="cluster", text_colour = "red") +
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))

# Slingshot --------------------------------------------------------------------

library(slingshot)
# TODO: Errors with MNN-corrected values
sce.sling <- slingshot(sce, reducedDim='PCA')
head(sce.sling$slingPseudotime_1)

# Setting up the colors.
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce.sling$slingPseudotime_1, breaks=100)]

# Creating a PCA plot.
plot(reducedDim(sce.sling, "PCA"), col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce.sling), lwd=2, col='black')

library(scater)
# sce.sling <- runUMAP(sce.sling, dimred="PCA")

# TODO: make ggcells robust to random crap in the colData().
# Also need to add a function to auto-generate a path.
# sce.sling$cell.type <- sce.sling$FACS <- NULL

library(viridis)
ggcells(sce.sling, mapping=aes(x=UMAP.1,
                               y=UMAP.2, col=slingPseudotime_1)) +
  geom_point() + scale_color_viridis()


sce.sling2 <- slingshot(sce, clusterLabels = sce$cluster, reducedDim='PCA')

plot(reducedDim(sce.sling2, "PCA"), col="grey80", pch=16, asp = 1)
lines(SlingshotDataSet(sce.sling2), lwd=2, col='black')

ggcells(sce.sling2, mapping=aes(x=UMAP.1,
                               y=UMAP.2, col=slingPseudotime_1)) +
  geom_point() + scale_color_viridis()
