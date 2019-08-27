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


#' plotFactorScatters2
#'
#' adapted from MOFA package: https://github.com/bioFAM/MOFA
#'
#' @param object MOFA object
#' @param factors
#' @param showMissing
#' @param color_by
#' @param name_color
#' @param shape_by
#' @param name_shape
#'
#' @return ggplot object
#' @export
#'
plotFactorScatters2 <- function(object, factors = "all", showMissing = TRUE, color_by = NULL,
                                name_color = "", shape_by = NULL, name_shape = "", col_values = NULL) {
  if (!is(object, "MOFAmodel")) {
    stop("'object' has to be an instance of MOFAmodel")
  }
  N <- object@Dimensions[["N"]]
  Z <- getFactors(object, factors = factors)
  factors <- colnames(Z)
  if (is.numeric(factors)) {
    factors <- factorNames(object)[factors]
  }
  else {
    if (paste0(factors, collapse = "") == "all") {
      factors <- factorNames(object)
    }
    else {
      stopifnot(all(factors %in% factorNames(object)))
    }
  }
  Z <- getFactors(object, factors = factors)
  tmp <- apply(Z, 2, var, na.rm = TRUE)
  if (any(tmp == 0)) {
    Z <- Z[, !tmp == 0]
    factors <- factors[!tmp == 0]
  }
  colorLegend <- TRUE
  if (!is.null(color_by)) {
    if (length(color_by) == 1 & is.character(color_by)) {
      if (name_color == "") {
        name_color <- color_by
      }
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData, rownames)
      if (color_by %in% Reduce(union, featureNames)) {
        viewidx <- which(vapply(featureNames, function(vnm) color_by %in%
            vnm, logical(1)))
        color_by <- TrainData[[viewidx]][color_by, ]
      }
      else if (is(object@InputData, "MultiAssayExperiment")) {
        color_by <- getCovariates(object, color_by)
      }
      else {
        stop("'color_by' was specified but it was not recognised, please read the documentation")
      }
    }
    else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
    }
    else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  }
  else {
    color_by <- rep(TRUE, N)
    colorLegend <- FALSE
  }
  shapeLegend <- TRUE
  if (!is.null(shape_by)) {
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if (name_shape == "") {
        name_shape <- shape_by
      }
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData, rownames)
      if (shape_by %in% Reduce(union, featureNames)) {
        viewidx <- which(vapply(featureNames, function(vnm) shape_by %in%
            vnm, logical(1)))
        shape_by <- TrainData[[viewidx]][shape_by, ]
      }
      else if (is(object@InputData, "MultiAssayExperiment")) {
        shape_by <- getCovariates(object, shape_by)
      }
      else {
        stop("'shape_by' was specified but it was not recognised, please read the documentation")
      }
    }
    else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    }
    else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
  }
  else {
    shape_by <- rep(TRUE, N)
    shapeLegend <- FALSE
  }
  if (!showMissing) {
    Z <- Z[!(is.na(color_by) | is.nan(color_by)), ]
    color_by <- color_by[!is.na(color_by)]
    shape_by <- shape_by[!is.na(shape_by)]
  }
  df <- as.data.frame(Z)
  colnames(df) <- paste0("LF", colnames(df))
  df <- cbind(df, color_by = color_by, shape_by = shape_by)
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if (length(unique(df$color_by)) < 5) {
    df$color_by <- as.factor(df$color_by)
  }
  main <- ""
  p <- ggplot(df, aes_string(
    x = colnames(df)[1], y = colnames(df)[2],
    color = "color_by", shape = "shape_by"
  )) + geom_point()

  p <- ggplot(df, aes_string(
    x = colnames(df)[1], y = colnames(df)[2],
    color = "color_by", shape = "shape_by"
  )) + geom_point()

  if (!is.null(col_values)) {
    p <- p + scale_color_manual("", values = col_values)
  }

  if (colorLegend | shapeLegend) {
    p <- p + theme(
      legend.title = element_text(
        size = 15,
        hjust = 0.5, color = "black"
      ), legend.position = "right",
      legend.direction = "vertical", legend.key = element_blank()
    )
    if (is.numeric(df$color_by)) {
      p <- p + scale_color_gradientn(colors = terrain.colors(10))
    }
    if (colorLegend) {
      p <- p + labs(color = name_color)
    }
    else {
      p <- p + guides(color = FALSE)
    }
    if (shapeLegend) {
      p <- p + labs(shape = name_shape)
    }
    else {
      p <- p + guides(shape = FALSE)
    }
    legend <- GGally::grab_legend(p)
  }
  else {
    legend <- NULL
  }
  p <- GGally::ggpairs(df,
    columns = colnames(df[, !colnames(df) %in% c("color_by", "shape_by")]),
    lower = list(continuous = "points"),
    diag = list(continuous = "barDiag"),
    upper = list(continuous = "points"),
    mapping = aes(color = color_by, shape = shape_by),
    title = main,
    legend = legend
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.direction = "vertical"
    )


  for (i in seq_len(p$nrow)) {
    for (j in seq_len(p$ncol)) {
      if (is.numeric(df$color_by)) {
        p[i, j] <- p[i, j] + scale_color_gradientn(colors = terrain.colors(10))
      } else if (!is.null(col_values)) {
        p[i, j] <- p[i, j] + scale_color_manual("", values = col_values)
        p[i, j] <- p[i, j] + scale_fill_manual("", values = col_values)
      }
    }
  }

  return(p)
}
