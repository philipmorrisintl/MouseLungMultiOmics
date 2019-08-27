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


GetOrtholog <-
function (Gene.Symbol, species_from = c("Hs", "Mm", "Rn")[1], 
    species_to = c("Hs", "Mm", "Rn")[2], source.homolog = c("HCOP", 
        "NCBI")[1], HomologTable = NULL) 
{
    if (is.null(HomologTable)) {
        HomologTable = GetHomologTable(source.homolog = source.homolog)
    }
    stopifnot(source.homolog %in% c("Bioconductor", "NCBI", "HCOP"))
    if (!paste0("from ", species_from, " to ", species_to) %in% 
        names(HomologTable)) 
        stop(species_from, " or ", species_to, " are not in names(HomologTable)")
    comment <- attr(HomologTable, "comment")
    HomologTable <- HomologTable[[paste0("from ", species_from, 
        " to ", species_to)]]
    res <- sapply(Gene.Symbol, function(x) {
        res.tmp = as.character(HomologTable[HomologTable$Gene.Symbol.x == 
            x, "Gene.Symbol.y"])
        if (length(res.tmp) == 1) {
            return(res.tmp)
        }
        else if (length(res.tmp) == 0) {
            return(NA)
        }
        else {
            if (toupper(x) %in% toupper(res.tmp)) {
                return(res.tmp[which(toupper(res.tmp) == toupper(x))][1])
            }
            else {
                return(res.tmp[1])
            }
        }
    })
    attr(res, "comment") <- comment
    return(res)
}
