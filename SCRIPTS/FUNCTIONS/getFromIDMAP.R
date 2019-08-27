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

getFromIDMAP <-
  function(idmap, genes, type = c(
               "all", "foldChange", "adj.pvalue",
               "p.value"
             )[1], fixed = TRUE, exact = FALSE, ...) {
    genes <- cleanGeneName(genes)
    idmap <- lapply(idmap, function(x) {
      x$nodeLabel <- cleanGeneName(x$nodeLabel)
      return(x)
    })
    if (exact == FALSE) {
      tmp0 <- lapply(idmap, function(x) {
        x[unlist(grep2(toupper(genes), cleanGeneName(toupper(x$nodeLabel)),
          fixed = fixed, ...
        )), , drop = FALSE]
      })
      if (type == "all") {
        return(tmp0)
      }
      if (type != "all") {
        tmp <- sapply(idmap, function(x) {
          x[unlist(grep2(toupper(genes), cleanGeneName(toupper(x$nodeLabel)),
            fixed = fixed
          )), colnames(x) == type]
        })
        if (class(tmp) != "matrix") {
          nm <- names(tmp)
          tmp <- matrix(unlist(tmp), nrow = 1)
          colnames(tmp) <- nm
        }
        rownames(tmp) <- tmp0[[1]]$nodeLabel
        return(tmp)
      }
    }
    else {
      tmp0 <- lapply(idmap, function(x) {
        x[toupper(x$nodeLabel) %in% cleanGeneName(toupper(genes))
          , ,
          drop = FALSE
        ]
      })
      if (type == "all") {
        return(tmp0)
      }
      if (type != "all") {
        tmp <- sapply(idmap, function(x) {
          x[
            toupper(x$nodeLabel) %in% cleanGeneName(toupper(genes)),
            colnames(x) == type
          ]
        })
        if (class(tmp) != "matrix") {
          nm <- names(tmp)
          tmp <- matrix(unlist(tmp), nrow = 1)
          colnames(tmp) <- nm
        }
        rownames(tmp) <- tmp0[[1]]$nodeLabel
        return(tmp)
      }
    }
  }
