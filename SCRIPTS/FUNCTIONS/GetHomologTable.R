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


#' This function GetHomologTable give the map table for homolog genes. Some core codes are from Bjoern Titz.
#' @export
#' @param build is a \code{character} \code{vector} with length 1. The latest version from NCBI website will be used if it is "current"; The latest saved version will be used if build is "saved"; Different build version will be used if build is "build*". for reproducibility, normal user can only use "saved" version.
#' @param species2code is an \code{list} which gives the taxonomy id for different species (http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/). If we add more species in species2code, then the mapping of the additional species will be added. For example, "Taxonomy ID" of zebrafish is 7955, so if we use species2code=list(Mm=10090, Rn = 10116, Hs = 9606, Zebrafish=7955), zebrafish will be added.
#' @param A saved version with file name useThisFile will be used if useThisFile is not NULL. Default is NULL.
#' @return A list.
#' @author Yang Xiang \email{yang.xiang@@pmi.com}; Bjoern Titz \email{bjorn.titz@@pmi.com}
#' @examples
#' # HomologTable <- GetHomologTable();
#' ## Or you can do for saved version
#' # data(HomologTable)
#' # if("from Hs to Mm" %in% names(HomologTable)) {stopifnot(c("C4b", "C4a") %in% as.vector(as.matrix(HomologTable$"from Hs to Mm"[HomologTable$"from Hs to Mm"$cluster_id == 36030 & HomologTable$"from Hs to Mm"$Gene.Symbol.x == "C4B", "Gene.Symbol.y"])))}
GetHomologTable <- function(build = c("current", "saved", "build67", "build68")[2],
                            species2code = list(Hs = 9606, Mm = 10090, Rn = 10116),
                            dir.save = getwd(),
                            useThisFile = NULL,
                            source.homolog = c("NCBI", "HCOP")[2]) {
  if (build != "saved") {
    warning("Only saved version is used for reproducibility of our analysis result. Please contact the developer or administrator for current (latest) or previous build.")
    build <- "saved" # To keep stable result, only the developer or administrator can use this option.
  }
  if (!is.null(useThisFile)) build <- "saved"
  ## Main section.

  dir.fun <- getwd()
  CO <- expand.grid(names(species2code), names(species2code))
  CO <- as.matrix(CO[CO[, 1] != CO[, 2], ])
  colnames(CO) <- c("from_species", "to_species")

  res.2 <- vector("list", nrow(CO))
  names(res.2) <- paste("from ", CO[, 1], " to ", CO[, 2], sep = "")

  if (build != "saved") {
    if (source.homolog == "HCOP") {
      cat("ERROR\n")
      cat("Use hcopCreateHomologTable() to create updated hcop table\n")
      stop("use different function to create hcop orthology table")
    }
    fn <- paste("ftp://ftp.ncbi.nih.gov/pub/HomoloGene/", build, "/homologene.data", sep = "")
    setwd(dir.save)
    fn.save <- paste(dir.save, "/homologene.data", sep = "")
    if (file.access(fn.save) == 0) unlink(fn.save)
    stopifnot(try(download.file(fn, fn.save, method = "wget")) != "try-error")
    fn <- fn.save

    res <- read.table(fn, sep = "\t", header = FALSE, as.is = TRUE, comment.char = "", quote = "")
    colnames(res) <- c("cluster_id", "species_id", "GeneID", "Gene.Symbol", "ncbi_prot_ref", "ncbi_seq")
    for (i in 1:nrow(CO)) { # i = c(1:nrow(CO))[3]
      species_from_id <- species2code[[CO[i, "from_species"]]]
      species_to_id <- species2code[[CO[i, "to_species"]]]
      tab_from <- res[res$species_id == species_from_id, c("cluster_id", "species_id", "GeneID", "Gene.Symbol")]
      tab_to <- res[res$species_id == species_to_id, c("cluster_id", "species_id", "GeneID", "Gene.Symbol")]

      tab_join <- base::merge(tab_from, tab_to, by = "cluster_id", all.x = TRUE)
      res.2[[i]] <- tab_join
    }

    HomologTable <- res.2

    ## save GeneIDTable

    save(HomologTable, file = paste(dir.save, "/HomologTable.",
      strsplit(as.character(Sys.time()), " ")[[1]][1], ".", build, ".rda",
      sep = ""
    ))
    cat("\nSave HomologTable to file ", paste(dir.save, "/HomologTable.",
      strsplit(as.character(Sys.time()), " ")[[1]][1], ".", build, ".rda",
      sep = ""
    ), "\n")
  } else { # useLatest == FALSE
    if (source.homolog == "NCBI") {
      HomologTable <- load2("../DATA/HomologTable.rda")
      return(HomologTable)
    } else if (source.homolog == "HCOP") {
      HCOP_ORTHOLOGY_TABLE <- load2("../DATA/HCOP_ORTHOLOGY_TABLE.rda")
      return(HCOP_ORTHOLOGY_TABLE)
    } else {
      stop("selected homology table not available")
    }
  }


  ## CP:
  if ("from Hs to Mm" %in% names(HomologTable)) {
    stopifnot(c("C4b", "C4a") %in% as.vector(as.matrix(HomologTable$"from Hs to Mm"[HomologTable$"from Hs to Mm"$cluster_id == 36030 & HomologTable$"from Hs to Mm"$Gene.Symbol.x == "C4B", "Gene.Symbol.y"])))
  }
  if ("from Hs to Mm" %in% names(HomologTable)) {
    stopifnot(c("LOC100910804", "Olr163", "LOC100910075", "Olr159") %in% HomologTable$"from Hs to Rn"[HomologTable$"from Hs to Rn"$cluster_id == "105161" & HomologTable$"from Hs to Rn"$Gene.Symbol.x == "OR52N4", "Gene.Symbol.y"])
  }

  ## Set back dir.
  setwd(dir.fun)

  ## Return
  return(HomologTable)
}
