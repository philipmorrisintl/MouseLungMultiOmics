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


pcaGG <-
function (X, group = NULL, cols = NULL, show = 1:2, main = "PCA Scores", 
    size = 1, size.text = 3, print = TRUE, ncol = NULL, panel = "Spectral", 
    annotate = FALSE) 
{
    if (is.null(group)) {
        group = factor(rep("", nrow(X)))
        cols = "navy"
    }
    if (is.null(cols)) {
        cols = colorRampPalette(brewer.pal(n = 7, panel))(nlevels(group))
    }
    X = as.matrix(X)
    if (is.null(rownames(X))) {
        rownames(X) = 1:nrow(X)
    }
    X = X[, complete.cases(t(X))]
    sv = svd(scale(X), nu = max(show), nv = max(show))
    scores = sv$u
    inertia = sv$d^2/sum(sv$d^2)
    colnames(scores) = paste0("PC", 1:ncol(scores), " (", round(100 * 
        inertia[1:ncol(scores)], 1), ")")
    varnm = paste0("PC", 1:ncol(scores), " (", round(100 * inertia[1:ncol(scores)], 
        1), "%)")
    D0 = data.frame(Id = rownames(X), Group = group, scores[, 
        show])
    colnames(D0)[-c(1:2)] = varnm[show]
    D0 = reshape2::melt(D0)
    D = NULL
    find_hull <- function(X) {
        ind = chull(X$x, X$y)
        return(data.frame(Id = rownames(X)[ind], X[ind, ]))
    }
    hulls = NULL
    for (i in 1:ncol(scores)) {
        for (j in 1:ncol(scores)) {
            if (i < j & all(c(i, j) %in% show)) {
                D2 = cbind(D0[D0$variable == varnm[i], ], D0[D0$variable == 
                  varnm[j], -c(1:2)])
                colnames(D2)[-c(1:2)] = c("Var1", "x", "Var2", 
                  "y")
                D = rbind(D, D2)
                tmp = cbind(D0[D0$variable == varnm[i], c("Group", 
                  "value")], D0[D0$variable == varnm[j], "value"])
                colnames(tmp) = c("Group", "x", "y")
                tmp = ddply(tmp, "Group", find_hull)
                tmp = data.frame(tmp[, c(1:2)], rep(varnm[i], 
                  nrow(tmp)), tmp[, "x"], rep(varnm[j], nrow(tmp)), 
                  tmp[, "y"])
                colnames(tmp) = colnames(D)
                hulls = rbind(hulls, tmp)
            }
        }
    }
    D$Vars = paste0(D$Var1, " vs. ", D$Var2)
    hulls$Vars = paste0(hulls$Var1, " vs. ", hulls$Var2)
    if (annotate == TRUE) {
        size = 0.001 * size
    }
    p = ggplot(data = D, aes(x = x, y = y, label = Id, fill = Group)) + 
        geom_point(size = size, shape = 21) + geom_hline(yintercept = 0, 
        colour = "darkgrey") + geom_vline(xintercept = 0, colour = "darkgrey") + 
        scale_colour_manual("", values = cols) + scale_fill_manual("", 
        values = cols) + labs(x = "", y = "") + ggtitle(main)
    if (annotate == TRUE) {
        p = p + geom_text(aes(colour = Group), size = size.text)
    }
    if (is.null(ncol)) {
        p = p + facet_wrap(~Vars)
    }
    else {
        p = p + facet_wrap(~Vars, ncol = ncol)
    }
    if (nlevels(group) == 1) {
        p = p + theme(legend.position = "none")
    }
    p = p + geom_polygon(data = hulls, aes(fill = Group), alpha = 0.2)
    if (print == TRUE) {
        print(p)
    }
    return(invisible(p))
}
