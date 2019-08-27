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


getIDMAPentry <-
  function(idmap, entry = "foldChange", GS = FALSE) {
    nr <- sapply(idmap, nrow)
    flag <- 0
    if (!all(nr == nr[1])) {
      flag <- 1
    }
    else {
      nm <- sapply(idmap, function(G) {
        G$nodeLabel
      })
      if (nr == 1) {
        nm <- matrix(nm, nrow = 1)
      }
      if (!all(apply(nm, 1, function(x) {
        all(x == x[1])
      }))) {
        flag <- 1
      }
    }
    if (flag == 1) {
      message("Re-arranging unequal sized idmap...")
      allids <- as.character(sort(unique(unlist(lapply(
        idmap,
        function(x) {
          x$nodeLabel
        }
      )))))
      idmap <- cleanIDMAP(idmap)
    }
    out <- sapply(idmap, function(x) x[, colnames(x) == entry])
    if (nrow(idmap[[1]]) == 1) {
      out <- matrix(out, nrow = 1)
      colnames(out) <- names(idmap)
    }
    rownames(out) <- idmap[[1]]$nodeLabel
    if (GS == TRUE) {
      rownames(out) <- cleanGeneName(as.character(idmap[[1]]$nodeLabel))
    }
    return(out)
  }
