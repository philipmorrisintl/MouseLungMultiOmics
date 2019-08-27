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


load2 <-
function (f, obj = NULL, verbose = TRUE) 
{
    env <- new.env()
    if (verbose == TRUE) {
        message("From File, ")
    }
    nm <- load(f, envir = env, verbose = verbose)
    if (verbose == TRUE) {
        message("\n")
    }
    res <- lapply(as.list(nm), function(n) {
        get(n, envir = env)
    })
    names(res) <- nm
    if (length(nm) == 1) {
        res <- res[[1]]
    }
    else {
        if (!is.null(obj)) {
            if (length(obj) == 1) {
                if (!obj %in% nm) {
                  stop("obj does not belong to loaded objects")
                }
                else {
                  res <- res[[which(names(res) == obj)]]
                }
            }
            else {
                stop("obj should be of length 1")
            }
        }
    }
    return(res)
}
