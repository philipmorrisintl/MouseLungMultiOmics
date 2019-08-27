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


getAttr <- function(x) {
  nm.attr <- names(attributes(x))
  nm.attr <- nm.attr[!nm.attr %in% c("names", "row.names")]
  attr0 <- NULL
  if (length(nm.attr) > 0) {
    attr0 <- vector("list", length(nm.attr))
    names(attr0) <- nm.attr
    for (n in 1:length(attr0)) {
      attr0[[n]] <- attr(x, nm.attr[n])
    }
    names(attr0) <- nm.attr
  }
  return(attr0)
}
