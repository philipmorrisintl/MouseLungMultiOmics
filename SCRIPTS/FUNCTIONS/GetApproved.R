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


#' This function GetApproved substitute synonyms with approved gene symbols. Note: it works for both character vector and idmap.
#' @export
#' @param Gene.Symbol is a \code{character} vector, or idmap.
#' @param species is a \code{string} which shouldd be "Hs", or "Mm", or "Rn".
#' @param nomatch is a \code{character} string.
#' @param ordered.Gene.Symbol is a \code{logical} string - returned result to be reordered.
#' @return A \code{character} vector if the input Gene.Symbol is a character vector, or and reordered idmap (by default) if Gene.Symbol is an idmap.
#' @author Yang Xiang \email{yang.xiang@@pmi.com}
GetApproved <- function(Gene.Symbol, species = c("Hs", "Mm", "Rn")[1], nomatch = c("Synonyms", "NA")[1], ordered.Gene.Symbol = TRUE) {
  useLatest <- FALSE

  if (is.list(Gene.Symbol)) { # Gene.Symbol is idmap.
    inputIsidmap <- TRUE
    idmap <- Gene.Symbol
    #if (!CheckNodeLabel(idmap)) stop("In function GetApproved, the first input is an idmap but with different nodeLabel for different contrasts. Please check its nodeLabel.")
    # Gene.Symbol <-
    Gene.Symbol <- gsub("^exp\\(|\\)$", "", as.character(idmap[[1]]$nodeLabel), ignore = TRUE)
  } else {
    inputIsidmap <- FALSE
  }

  geneIDTable <- GetGeneIDTable(species = species, useLatest = useLatest)
  approvedTable <- GetApprovedTable(species = species, useLatest = useLatest)
  approvedSymbol <- geneIDTable[geneIDTable[, "type"] == "Approved", "Gene.Symbol"]
  ind.syn <- which(!Gene.Symbol %in% approvedSymbol)
  if (length(ind.syn) != 0) {
    syn <- Gene.Symbol[ind.syn]
    syn.mem <- syn
    ind <- match(syn, approvedTable[, "Synonyms"])
    if (nomatch == "NA") {
      syn <- approvedTable[, "Approved.Symbol"][ind]
    } else {
      syn[!is.na(ind)] <- approvedTable[, "Approved.Symbol"][ind[!is.na(ind)]]
    }
    Gene.Symbol[ind.syn] <- syn
  } else {
    cat("All the gene symbols are approved.\n")
  }

  if (inputIsidmap) {
    idmap.new <- lapply(idmap, function(x) {
      x$nodeLabel <- paste("exp(", Gene.Symbol, ")", sep = "")
      if (ordered.Gene.Symbol) x <- x[order(x$nodeLabel), ]
      x$nodeLabel <- factor(x$nodeLabel)
      return(x)
    })
    return(idmap.new)
  }

  return(Gene.Symbol)
}
