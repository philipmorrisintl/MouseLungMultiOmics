# Copyright 2019 Philip Morris Products SA
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


getIDMAPheatmapGG <-
  function(idmap, genes, gene.categories = NULL, species = NULL,
             BW = FALSE, pvthresh = c(0.05, 0.01), pvthresh.symbols = c(
               4,
               8
             ), symbol_size = 2, pvtype = c("p.value", "adj.p.value")[2],
             subtitle = NULL, key.title = "log2(FC)", title = "log2(Fold-change) Heatmap",
             ...) {
    useSymbols <- FALSE
    if (is.numeric(pvthresh.symbols)) {
      useSymbols <- TRUE
    }
    if (length(pvthresh) != length(pvthresh.symbols)) {
      pvsymb <- sapply(1:length(pvthresh), function(i) {
        paste(rep(pvthresh.symbols[1], i), collapse = "")
      })
    }
    else {
      pvsymb <- pvthresh.symbols[order(pvthresh, decreasing = TRUE)]
    }
    pvtype2 <- c("p-value", "fdr")[match(pvtype, c(
      "p.value",
      "adj.p.value"
    ))]
    message(pvthresh[order(pvthresh, decreasing = TRUE)])
    message(pvsymb)
    fcs <- getIDMAPentry(idmap)
    pvs <- getIDMAPentry(idmap, entry = pvtype)
    gin0 <- genes
    gin <- gin0[gin0 %in% rownames(fcs)]
    if (is.null(gene.categories)) {
      gene.categories <- rep("", length(gin))
    }
    else {
      gene.categories <- factor(gene.categories[gin0 %in% gin])
    }
    gr <- factor(gene.categories)
    fcs <- fcs[match(gin, rownames(fcs)), , drop = FALSE]
    pvs <- pvs[match(gin, rownames(pvs)), , drop = FALSE]
    if (max(table(paste(rownames(fcs), gr))) > 1) {
      t0 <- table(paste(rownames(fcs), gr))
      gnm <- rownames(fcs)
      fcs <- fcs[!duplicated(paste(gnm, gr)), , drop = FALSE]
      pvs <- pvs[!duplicated(paste(gnm, gr)), , drop = FALSE]
      gr <- gr[!duplicated(paste(gnm, gr))]
      warning(paste0(
        "Some rownames are not unique (within group, if any):",
        paste(names(t0)[t0 > 1], collapse = ", ")
      ))
    }
    if (BW == FALSE) {
      txt <- matrix(NA, nrow(pvs), ncol(pvs))
      pvthresh <- pvthresh[order(pvthresh, decreasing = TRUE)]
      for (k in 1:length(pvthresh)) {
        txt[pvs < pvthresh[k]] <- pvsymb[k]
      }
      symlab <- paste0("<", pvthresh)
      names(symlab) <- pvsymb
      ImagePlotGG(fcs,
        textmat = txt, useSymbols = useSymbols,
        symbol_labels = symlab, symbol_key_title = pvtype2,
        group = gr, subtitle = subtitle, symbol_size = symbol_size,
        key.title = key.title, title = title, ...
      )
    }
    else {
      if (length(pvthresh) > 1) {
        warning("In BW only two thresholds are allowed.")
        pvthresh <- pvthresh[1]
      }
      pvsymb <- c(24, 25)
      txt <- matrix(NA, nrow(pvs), ncol(pvs))
      txt[pvs < pvthresh[1] & fcs <= 0] <- pvsymb[2]
      txt[pvs < pvthresh[1] & fcs > 0] <- pvsymb[1]
      symlab <- c(
        paste0("log2(fc)>0 & ", pvtype2, " <", pvthresh),
        paste0("log2(fc)<0 & ", pvtype2, " <", pvthresh)
      )
      names(symlab) <- pvsymb
      sub0 <- subtitle
      ImagePlotGG(abs(fcs),
        textmat = txt, useSymbols = useSymbols,
        symbol_labels = symlab, symbol_key_title = pvtype2,
        group = gr, subtitle = subtitle, symbol_size = symbol_size,
        key.title = paste0("|", key.title, "|"), title = title,
        col.scale = c(
          "gray90", "gray70", "gray50", "grey40",
          "gray30"
        ), ...
      )
    }
  }
