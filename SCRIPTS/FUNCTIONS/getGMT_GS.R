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

getGMT_GS <- function(gmt_file_name) {
    gmt <- readr::read_delim(gmt_file_name, delim = "\t", col_names = FALSE) %>% 
        as.data.frame()
    gsc <- list()
    for (i in seq_len(nrow(gmt))) {
        genes <- unlist(gmt[i,3:ncol(gmt), drop = TRUE])
        genes <- genes[!is.na(genes)]
        gs_name <- gmt[i,1]
        gsc[[gs_name]] <- unname(genes)
    }
    return(gsc)
}