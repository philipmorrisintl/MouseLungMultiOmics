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


getLayoutMatrix <-
function (which_order_volcano, norow, nocol) 
{
    layout = matrix(NA, nrow = norow, ncol = nocol)
    count = 1
    for (i in 1:norow) {
        for (k in 1:nocol) {
            if (count > length(which_order_volcano)) {
                pos = count
            }
            else {
                pos = which_order_volcano[count]
            }
            if (pos == -1) {
                pos = max(which_order_volcano, layout, na.rm = TRUE) + 
                  1
            }
            layout[i, k] = pos
            count = count + 1
        }
    }
    return(layout)
}
