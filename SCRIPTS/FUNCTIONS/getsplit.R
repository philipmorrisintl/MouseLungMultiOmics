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


getsplit <-
function (nm0, splitarg, k = 1, remove = FALSE, fixed = TRUE, 
    last = FALSE, last.n = 1) 
{
    nm0 <- as.character(nm0)
    if (last == FALSE) {
        if (remove == FALSE) {
            y <- sapply(strsplit(nm0, splitarg, fixed = fixed), 
                function(x) {
                  paste(x[k], collapse = splitarg, sep = "")
                })
        }
        if (remove == TRUE) {
            y <- sapply(strsplit(nm0, splitarg, fixed = fixed), 
                function(x) {
                  paste(x[-k], collapse = splitarg, sep = "")
                })
        }
    }
    else {
        if (remove == FALSE) {
            y <- sapply(strsplit(nm0, splitarg, fixed = fixed), 
                function(x) {
                  paste(x[(length(x) - last.n + 1):length(x)], 
                    collapse = splitarg, sep = "")
                })
        }
        if (remove == TRUE) {
            y <- sapply(strsplit(nm0, splitarg, fixed = fixed), 
                function(x) {
                  paste(x[-c((length(x) - last.n + 1):length(x))], 
                    collapse = splitarg, sep = "")
                })
        }
    }
    names(y) <- names(nm0)
    return(y)
}
