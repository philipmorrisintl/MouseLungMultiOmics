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


getGS <-
  function(name = NULL) {
    fp <- "../DATA/"
    if (is.null(name)) {
      cat("Possible (but not all ) values for name are:", fill = TRUE)
      nms <- gsub(".rda", "", dir(fp, pattern = "rda$"), fixed = TRUE)
      nms <- nms[!nms %in% c(
        "genesets.as.hyps.hs",
        "MarkerPanels", "GSAPM.example"
      )]
      cat(paste(nms, collapse = "\n"), fill = TRUE)
      return(nms)
    }
    else {
      GS <- load2(paste0(fp, "/", name, ".rda"))
      if (is.data.frame(GS[[1]])) {
        GS <- lapply(GS, function(x) unique(as.character(x$nodeLabel)))
      }
      return(GS)
    }
  }
