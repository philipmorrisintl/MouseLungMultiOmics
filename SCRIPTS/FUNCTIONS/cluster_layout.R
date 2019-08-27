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


cluster_layout <- function(graph, sel_weights = c(0.2, 1), seed = 98234) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  g2 <- graph

  cl_member <- V(graph)$group
  names(cl_member) <- V(graph)$name

  e_mat <- get.edges(g2, E(g2))
  e_mat[, 1] <- V(g2)$name[as.numeric(e_mat[, 1])]
  e_mat[, 2] <- V(g2)$name[as.numeric(e_mat[, 2])]
  e_mat <- as.data.frame(e_mat)
  e_mat$C1 <- cl_member[e_mat[, 1]]
  e_mat$C2 <- cl_member[e_mat[, 2]]
  eweights <- ifelse(e_mat$C1 == e_mat$C2, sel_weights[1], sel_weights[2])
  glayout2 <- layout.fruchterman.reingold(g2, weights = eweights)

  glayout <- glayout2[match(V(graph)$name, V(g2)$name), ]

  return(glayout)
}
