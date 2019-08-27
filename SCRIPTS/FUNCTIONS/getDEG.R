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


getDEG <-
  function(dL, fdr = 0.05, fcth = 0, full = TRUE, pvtype = c(
               "fdr",
               "raw"
             )[1]) {
    if (pvtype == "fdr") {
      if (full == FALSE) {
        return(lapply(dL, function(x) {
          as.character(x$nodeLabel[!is.na(x$adj.p.value) &
            x$adj.p.value < fdr & abs(x$foldChange) > fcth])
        }))
      }
      if (full == TRUE) {
        return(lapply(dL, function(x) {
          x[!is.na(x$adj.p.value) & x$adj.p.value < fdr &
            abs(x$foldChange) > fcth, , ]
        }))
      }
    }
    if (pvtype == "raw") {
      if (full == FALSE) {
        return(lapply(dL, function(x) {
          as.character(x$nodeLabel[!is.na(x$adj.p.value) &
            x$p.value < fdr & abs(x$foldChange) > fcth])
        }))
      }
      if (full == TRUE) {
        return(lapply(dL, function(x) {
          x[!is.na(x$adj.p.value) & x$p.value < fdr & abs(x$foldChange) >
            fcth, , ]
        }))
      }
    }
  }
