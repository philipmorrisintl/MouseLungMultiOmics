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


barchartNonOmicsMouse <- function(dat_comb, ep, laby, cols = cols_mouse) {
  dat <- dat_comb[dat_comb$EndpointMapped == ep, ]
  endpoint_name <- getsplit(ep, "|", 1:3)

  dat$Group <- factor(dat$Group, levels = levels(dat$Group))
  dat$y_range <- max(dat$Mean, na.rm = TRUE) - min(dat$Mean, na.rm = TRUE)

  p <- ggplot(aes(x = Group, y = Mean, fill = Group, label = p.value_label), data = dat)
  p <- p + geom_bar(stat = "identity", width = 0.6)
  p <- p + scale_fill_manual("", values = cols, guide = FALSE)
  p <- p + geom_errorbar(aes(ymin = Mean - StdErr, ymax = Mean + StdErr), width = .2)
  p <- p + geom_text(aes(y = Mean + y_range / 20, x = Group), colour = "red", size = 11)
  # p <- p + labs(x = "", y = laby, title = endpoint_name)
  p <- p + labs(x = "", y = laby, title = "")
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 12))
  p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, colour = "black", size = 12))
  p <- p + theme(plot.title = element_text(size = 12, face = "bold"))
  p <- p + theme(axis.title = element_text(size = 12, face = "bold"))
  p <- p + theme(strip.background = element_blank())
  p <- p + theme(strip.text = element_text(colour = "black", size = 12, face = "bold"))
  p <- p + geom_hline(yintercept = 0, colour = "black")
  p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))
  return(p)
}
