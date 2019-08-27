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


#' Clean an IDMAPS object
#'
#' Ensure that contrasts are consistent each others in terms
#' of genes symbols, ordering...
#'
#' @export
#' @import BiocGenerics
#' @param dL A \code{list} correponding to an IDMAPS R object.
#' @param upper \code{boolean}. TRUE to convert all 'nodeLabel' to uppercase, FAlSE by default.
#' @param replace.na \code{boolean}. TRUE to replace NAs by a numeric value.
#' @param replacement.value \code{numeric}. The numeric value to replace NAs.
#' @return A \code{list} correponding to an IDMAPS R object cleaned.
#' @author Florian Martin
cleanIDMAP <- function(dL, upper = FALSE, replace.na = TRUE, replacement.value = 0) {
  attr0 <- getAttr(dL)
  attr0.idmap <- lapply(dL, getAttr)
  in.col <- Reduce(intersect, lapply(dL, colnames))
  if (!all(c("foldChange", "t", "p.value", "adj.p.value") %in% in.col)) {
    stop("colnames of idmap entries not similar enough!")
  }
  dL <- lapply(dL, function(X) return(X[, colnames(X) %in% in.col]))

  if (upper == TRUE) {
    dL <- lapply(dL, function(x) {
      y <- x
      y$nodeLabel <- factor(toupper(y$nodeLabel))
      return(y)
    })
  } else {
    dL <- lapply(dL, function(x) {
      y <- x
      y$nodeLabel <- factor(y$nodeLabel)
      return(y)
    })
  }
  dL <- lapply(dL, function(x) {
    y <- x[!duplicated(as.character(x$nodeLabel)), ]
    return(y)
  })
  dL <- lapply(dL, function(x) {
    y <- x[order(as.character(x$nodeLabel)), ]
    return(y)
  })
  nm <- lapply(dL, function(G) {
    G$nodeLabel
  })
  ok <- sapply(nm, function(n) {
    ok <- TRUE
    if (length(n) != length(nm[[1]])) {
      ok <- FALSE
    } else {
      ok <- all(sort(as.character(n)) == sort(as.character(nm[[1]])))
    }
    return(ok)
  })
  if (!all(ok == TRUE)) {
    warning("Issue with ids in dL. Cleaning process is performed")
    allids <- as.character(sort(unique(unlist(lapply(dL, function(x) {
      as.character(x$nodeLabel)
    })))))
    if (upper) {
      allids <- toupper(allids[!duplicated(toupper(allids))])
    }
    dL <- lapply(dL, function(x) {
      y <- x
      if (!"p.value" %in% colnames(y)) {
        pv0 <- 0 * y$foldChange
        s <- y$foldChange / y$t
        s[y$t == 0] <- 0
        y$p.value <- 2 * pt(abs(y$foldChange) / (s + 1e-12), y$df, lower.tail = FALSE)
      }
      matNA <- matrix(NA, ncol = ncol(y), nrow = length(allids[!allids %in%
        x$nodeLabel]))
      matNA[, colnames(y) == "nodeLabel"] <- allids[!allids %in% x$nodeLabel]
      colnames(matNA) <- colnames(y)
      y <- as.data.frame(rbind(y, as.data.frame(matNA)))
      y <- y[order(as.character(y$nodeLabel)), ]
      if (replace.na == TRUE & is.na(replacement.value)) {
        stop("replacement.value is NA and relpace.na is TRUE!")
      }
      if (!is.na(replacement.value) & replace.na == TRUE) {
        y$foldChange[is.na(y$foldChange)] <- replacement.value
        y$p.value[is.na(y$p.value)] <- 1 - replacement.value
        y$adj.p.value[is.na(y$adj.p.value)] <- 1 - replacement.value
        y$t[is.na(y$t)] <- replacement.value
        if ("df" %in% colnames(y)) {
          y$df[is.na(y$df)] <- replacement.value
        }
      }
      y[, colnames(y) %in% c("foldChange", "t", "df", "p.value", "adj.p.value")] <- apply(
        y[
          ,
          colnames(y) %in% c("foldChange", "t", "df", "p.value", "adj.p.value")
        ],
        2, as.numeric
      )
      return(y)
    })

    nm <- sapply(dL, function(G) {
      G$nodeLabel
    })
    if (!all(apply(nm, 1, function(x) all(x == x[1])))) {
      stop("Issue with ids in dL")
    }
  }
  for (k in 1:length(dL)) {
    rownames(dL[[k]]) <- dL[[k]]$nodeLabel
    dL[[k]] <- setAttr(dL[[k]], attr0.idmap[[k]])
  }
  if (length(attr0) > 0) {
    dL <- setAttr(dL, attr0)
  }
  return(dL)
}
