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

getLimmaResults <-
function (X, I = NULL, contrast = NULL, covariates = NULL, model.formula = NULL, 
    model.type = c("multiple", "unique")[1], ps2gs = FALSE, gsmap = c("HG-U133_Plus_2", 
        "HG-U133A", "HG-U133B", "Rat230_2", "Mouse430_2", "HT_MG-430_PM", 
        "GPL6883", "GPL2986")[1], padj.method = "fdr", data.type = "affy", 
    collapsing.method = "contrast_data", verbose = TRUE, fn.pset2gs = c("pset2gs", 
        "pset2gs_new")[2], include.ordinary.pvalues = FALSE, 
    fitOnly = FALSE, keepLimmaFit = TRUE, ...) 
{
    match.arg(model.type, c("multiple", "unique"))
    species <- "Hs"
    if (gsmap %in% c("Mouse430_2", "HT_MG-430_PM")) {
        species = "Mm"
    }
    if (gsmap %in% c("Rat230_2")) {
        species = "Rn"
    }
    match.arg(species, c("Hs", "Mm", "Rn"))
    if (!is.null(contrast)) {
        if (is.data.frame(contrast)) {
            contrast = as.matrix(contrast)
        }
        if (!is.matrix(contrast)) {
            stop("contrast is not NULL and is neither a data.frame nor a matrix")
        }
    }
    keepFit = vector("list")
    if (!is.null(model.formula)) {
        if (data.type == "count") {
            if (verbose == TRUE) {
                cat("Beware EXPERIMENTAL!", fill = TRUE)
            }
            nf = calcNormFactors(X)
            X = voom(X, mm0, plot = TRUE, lib.size = colSums(X) * 
                nf)
            X$genes = rownames(X)
        }
        lm.fit1 = limma::lmFit(X, model.matrix(as.formula(model.formula), 
            data = data.frame(covariates)), ...)
        if (!is.null(contrast)) {
            if (verbose == TRUE) {
                cat("Contrast matrix used in the context of the formula.", 
                  fill = TRUE)
            }
            if (!all(rownames(contrast) %in% colnames(lm.fit1$coefficients))) {
                stop(paste("Rownames of the contrast matrix should be in ", 
                  paste(colnames(lm.fit1$coefficients), collapse = ",")))
            }
            lm.fit1 <- contrasts.fit(lm.fit1, contrasts = contrast)
        }
        if (fitOnly == TRUE) {
            fit.eb1 = lm.fit1
            nmAdd = c("df.prior", "s2.prior", "var.prior", "proportion", 
                "s2.post", "t", "df.total", "p.value", "lods", 
                "F", "F.p.value")
            for (nam in nmAdd) {
                fit.eb1[[nam]] = matrix(NA, nrow(lmfit$coefficients), 
                  ncol(lmfit$coefficients))
                rownames(fit.eb1[[nam]]) = rownames(lmfit$coefficients)
                colnames(fit.eb1[[nam]]) = colnames(lmfit$coefficients)
            }
        }
        else {
            fit.eb1 = eBayes(lm.fit1)
        }
        if (keepLimmaFit == TRUE) {
            fit.eb1$data = X
            keepFit = fit.eb1
        }
        res = vector("list", ncol(fit.eb1$p.value))
        names(res) = colnames(fit.eb1$p.value)
        if (verbose == TRUE) {
            cat(colnames(fit.eb1$coefficients), fill = TRUE)
        }
        for (k in 1:length(res)) {
            if (!include.ordinary.pvalues) {
                res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                  foldChange = fit.eb1$coefficients[, k], adj.p.value = p.adjust(fit.eb1$p.value[, 
                    k], method = padj.method), p.value = fit.eb1$p.value[, 
                    k], t = fit.eb1$t[, k], df = fit.eb1$df.prior + 
                    fit.eb1$df.residual, A = rowMeans(X), TTT = rowMeans(X) * 
                    0 - 1, CTRL = rowMeans(X) * 0 - 1)
            }
            else {
                ord.t = fit.eb1$coef[, k]/fit.eb1$stdev.unscaled[, 
                  k]/fit.eb1$sigma
                ord.p = 2 * pt(abs(ord.t), df = fit.eb1$df.residual, 
                  lower.tail = FALSE)
                ord.adjp = p.adjust(ord.p, method = padj.method)
                res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                  foldChange = fit.eb1$coefficients[, k], adj.p.value = p.adjust(fit.eb1$p.value[, 
                    k], method = padj.method), p.value = fit.eb1$p.value[, 
                    k], t = fit.eb1$t[, k], df = fit.eb1$df.prior + 
                    fit.eb1$df.residual, A = rowMeans(X, na.rm = TRUE), 
                  TTT = rowMeans(X, na.rm = TRUE) * 0 - 1, CTRL = rowMeans(X, 
                    na.rm = TRUE) * 0 - 1, ord.t = ord.t, ord.p.value = ord.p, 
                  ord.adj.p.value = ord.adjp)
            }
        }
    }
    if (is.null(model.formula)) {
        I = factor(as.character(I))
        if (!all(as.vector(as.matrix(contrast)) %in% I)) {
            stop("I must contain all the constrasts involved in contrast argument.")
        }
        if (verbose == TRUE) {
            cat("Computing Limma Models...")
        }
        N = nrow(contrast)
        flagN = 0
        res = vector("list", N)
        names(res) = contrast[1:N, 1]
        if (model.type == "unique") {
            mm = mm0 = model.matrix(~I - 1, data.frame(I))
            if (!is.null(covariates)) {
                if (!is.data.frame(covariates)) {
                  stop("covariates must be a dataframe")
                }
                colnames(covariates) = paste("ZZ_", colnames(covariates), 
                  sep = "")
                mm0 = model.matrix(~1 + ., data = data.frame(I, 
                  covariates))
            }
            if (data.type == "count") {
                nf = calcNormFactors(X)
                X = voom(X, mm0, plot = TRUE, lib.size = colSums(X) * 
                  nf)
                X$genes = rownames(X)
            }
            lm.fit1 <- limma::lmFit(X, mm0)
            contrast2 = matrix(0, ncol = nrow(contrast), nrow = ncol(mm))
            rownames(contrast2) = colnames(mm)
            colnames(contrast2) = paste("CONT", contrast[, 1], 
                sep = "")
            for (k in 1:nrow(contrast)) {
                contrast2[rownames(contrast2) == paste("I", contrast[k, 
                  2], sep = ""), colnames(contrast2) == paste("CONT", 
                  contrast[k, 1], sep = "")] = -1
                contrast2[rownames(contrast2) == paste("I", contrast[k, 
                  1], sep = ""), colnames(contrast2) == paste("CONT", 
                  contrast[k, 1], sep = "")] = 1
            }
            if (!is.null(covariates)) {
                tmp = matrix(0, ncol = ncol(contrast2), nrow = length(which(!colnames(lm.fit1$coefficients) %in% 
                  rownames(contrast2))))
                rownames(tmp) = colnames(lm.fit1$coefficients)[!colnames(lm.fit1$coefficients) %in% 
                  rownames(contrast2)]
                colnames(tmp) = colnames(contrast2)
                contrast2 = rbind(contrast2, tmp)
            }
            contrast2 = contrast2[match(colnames(lm.fit1$coefficients), 
                rownames(contrast2)), , drop = FALSE]
            lmfit <- contrasts.fit(lm.fit1, contrasts = contrast2)
            if (fitOnly == TRUE) {
                fit.eb1 = lmfit
                nmAdd = c("df.prior", "s2.prior", "var.prior", 
                  "proportion", "s2.post", "t", "df.total", "p.value", 
                  "lods", "F", "F.p.value")
                tmpAdd = matrix(NA, nrow(lmfit$coefficients), 
                  ncol(lmfit$coefficients))
                rownames(tmpAdd) = rownames(lmfit$coefficients)
                colnames(tmpAdd) = colnames(lmfit$coefficients)
                for (nam in nmAdd) {
                  fit.eb1[[nam]] = tmpAdd
                }
            }
            else {
                fit.eb1 <- eBayes(lmfit)
            }
            if (keepLimmaFit == TRUE) {
                fit.eb1$data = X
                keepFit = fit.eb1
            }
            if (verbose == TRUE) {
                cat(colnames(fit.eb1$p.value), fill = TRUE)
            }
            for (k in 1:N) {
                ttt = c(1:length(I))[I == contrast[k, 1]]
                ctrl = c(1:length(I))[I == contrast[k, 2]]
                if (length(ttt) == 1) {
                  Xttt = matrix(X[, ttt], ncol = 1)
                }
                else {
                  Xttt = X[, ttt]
                }
                if (length(ctrl) == 1) {
                  Xctrl = matrix(X[, ctrl], ncol = 1)
                }
                else {
                  Xctrl = X[, ctrl]
                }
                if (!include.ordinary.pvalues) {
                  res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                    foldChange = fit.eb1$coefficients[, k], adj.p.value = p.adjust(fit.eb1$p.value[, 
                      k], method = padj.method), p.value = fit.eb1$p.value[, 
                      k], t = fit.eb1$t[, k], df = fit.eb1$df.prior + 
                      fit.eb1$df.residual, A = rowMeans(X, na.rm = TRUE), 
                    TTT = rowMeans(Xttt, na.rm = TRUE), CTRL = rowMeans(Xctrl, 
                      na.rm = TRUE))
                }
                else {
                  ord.t = fit.eb1$coef[, k]/fit.eb1$stdev.unscaled[, 
                    k]/fit.eb1$sigma
                  ord.p = 2 * pt(abs(ord.t), df = fit.eb1$df.residual, 
                    lower.tail = FALSE)
                  ord.adjp = p.adjust(ord.p, method = padj.method)
                  res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                    foldChange = fit.eb1$coefficients[, k], adj.p.value = p.adjust(fit.eb1$p.value[, 
                      k], method = padj.method), p.value = fit.eb1$p.value[, 
                      k], t = fit.eb1$t[, k], df = fit.eb1$df.prior + 
                      fit.eb1$df.residual, A = rowMeans(X, na.rm = TRUE), 
                    TTT = rowMeans(Xttt, na.rm = TRUE), CTRL = rowMeans(Xctrl, 
                      na.rm = TRUE), ord.t = ord.t, ord.p.value = ord.p, 
                    ord.adj.p.value = ord.adjp)
                }
            }
        }
        else {
            for (k in 1:N) {
                if (verbose == TRUE) {
                  cat(paste("\n###### Estimating", contrast[k, 
                    1], "vs", contrast[k, 2], "...\n"))
                }
                ttt = c(1:length(I))[I == contrast[k, 1]]
                ctrl = c(1:length(I))[I == contrast[k, 2]]
                fact = c(paste("B__", I[ttt], sep = ""), paste("A__", 
                  I[ctrl], sep = ""))
                Xtmp = X[, c(ttt, ctrl)]
                if (is.null(covariates)) {
                  design.mat1 <- model.matrix(~1 + TTT, data = data.frame(TTT = fact))
                }
                if (!is.null(covariates)) {
                  if (!is.data.frame(covariates)) {
                    stop("covariates must be a dataframe")
                  }
                  if (k == 1) {
                    colnames(covariates) = paste("ZZ_", colnames(covariates), 
                      sep = "")
                  }
                  cov = covariates[c(ttt, ctrl), , drop = FALSE]
                  for (h in 1:ncol(cov)) {
                    if (class(cov[, h]) == "factor") {
                      cov[, h] = factor(cov[, h])
                    }
                  }
                  colnames(cov) = colnames(covariates)
                  design.mat1 = model.matrix(as.formula(paste("~1", 
                    paste(colnames(cov), collapse = "+"), "ATTT", 
                    sep = "+")), data = data.frame(ATTT = fact, 
                    cov))
                  if (verbose == TRUE) {
                    cat(paste("~1", paste(colnames(cov), collapse = "+"), 
                      "ATTT", sep = "+"), fill = TRUE)
                  }
                }
                if (data.type == "count") {
                  nf = calcNormFactors(Xtmp)
                  Xtmp = voom(Xtmp, design.mat1, plot = TRUE)
                  Xtmp$genes = rownames(X)
                }
                lm.fit1 <- try(limma::lmFit(Xtmp, design.mat1), 
                  silent = TRUE)
                if (fitOnly == TRUE) {
                  fit.eb1 = lm.fit1
                  nmAdd = c("df.prior", "s2.prior", "var.prior", 
                    "proportion", "s2.post", "t", "df.total", 
                    "p.value", "lods", "F", "F.p.value")
                  tmpAdd = matrix(NA, nrow(lm.fit1$coefficients), 
                    ncol(lm.fit1$coefficients))
                  rownames(tmpAdd) = rownames(lm.fit1$coefficients)
                  colnames(tmpAdd) = colnames(lm.fit1$coefficients)
                  for (nam in nmAdd) {
                    fit.eb1[[nam]] = tmpAdd
                  }
                }
                else {
                  fit.eb1 <- try(eBayes(lm.fit1), silent = TRUE)
                }
                if (keepLimmaFit == TRUE) {
                  fit.eb1$data = Xtmp
                  keepFit[[length(keepFit) + 1]] = fit.eb1
                  names(keepFit)[length(keepFit)] = paste(contrast[k, 
                    ], collapse = " vs ")
                }
                if (class(lm.fit1) == "try-error" || class(fit.eb1) == 
                  "try-error") {
                  warning("lmFit and/or eBayes fitting fail for ", 
                    contrast[k, 1])
                  next
                }
                if (verbose == TRUE) {
                  cat(colnames(fit.eb1$coefficients), fill = TRUE)
                }
                if (data.type == "count") {
                  Xtmp = Xtmp$E
                }
                if (length(ttt) == 1) {
                  Xttt = matrix(X[, ttt], ncol = 1)
                }
                else {
                  Xttt = X[, ttt]
                }
                if (length(ctrl) == 1) {
                  Xctrl = matrix(X[, ctrl], ncol = 1)
                }
                else {
                  Xctrl = X[, ctrl]
                }
                nc = grep("TTTB__", colnames(fit.eb1$coefficients))
                if (verbose == TRUE) {
                  cat(colnames(fit.eb1$coefficients)[nc], fill = TRUE)
                }
                if (!include.ordinary.pvalues) {
                  res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                    foldChange = fit.eb1$coefficients[, nc], 
                    adj.p.value = p.adjust(fit.eb1$p.value[, 
                      nc], method = padj.method), p.value = fit.eb1$p.value[, 
                      nc], t = fit.eb1$t[, nc], df = fit.eb1$df.prior + 
                      fit.eb1$df.residual, A = rowMeans(Xtmp, 
                      na.rm = TRUE), TTT = rowMeans(Xttt, na.rm = TRUE), 
                    CTRL = rowMeans(Xctrl, na.rm = TRUE))
                }
                else {
                  ord.t = fit.eb1$coef[, nc]/fit.eb1$stdev.unscaled[, 
                    nc]/fit.eb1$sigma
                  ord.p = 2 * pt(abs(ord.t), df = fit.eb1$df.residual, 
                    lower.tail = FALSE)
                  ord.adjp = p.adjust(ord.p, method = padj.method)
                  res[[k]] = data.frame(nodeLabel = rownames(fit.eb1$coefficients), 
                    foldChange = fit.eb1$coefficients[, nc], 
                    adj.p.value = p.adjust(fit.eb1$p.value[, 
                      nc], method = padj.method), p.value = fit.eb1$p.value[, 
                      nc], t = fit.eb1$t[, nc], df = fit.eb1$df.prior + 
                      fit.eb1$df.residual, A = rowMeans(Xtmp), 
                    TTT = rowMeans(Xttt, na.rm = TRUE), CTRL = rowMeans(Xctrl, 
                      na.rm = TRUE), ord.t = ord.t, ord.p.value = ord.p, 
                    ord.adj.p.value = ord.adjp)
                }
            }
        }
    }
    if (verbose == TRUE) {
        cat("\n")
    }
    ind <- which(!sapply(res, is.null))
    if (!all(c(1:length(res)) %in% ind)) {
        warning("\ncontrasts are removed: ", paste(names(res)[!c(1:length(res)) %in% 
            ind], sep = " ", collapse = ", "), "\n")
        res <- res[ind]
    }
    res0 = res
    if (ps2gs == TRUE) {
        if (verbose == TRUE) {
            cat("Converting results to Gene Symbol...")
            cat(paste("map:", gsmap, "\n"))
        }
        nmres = names(res)
        nmcolres = names(res[[1]])
        adjust.p.values = FALSE
        statistics.type = "other"
        if ("t" %in% nmcolres) {
            adjust.p.values = TRUE
            statistics.type = "t"
        }
        res.load <- do.call("library", list("RConferoMapping"))
        if (fn.pset2gs == "pset2gs") {
            res = pset2gs(res, map = gsmap, species = species, 
                collapsing.method = collapsing.method, adjust.p.values = adjust.p.values, 
                statistics.type = statistics.type)
        }
        else {
            res = pset2gs_new(res, map = gsmap, species = species, 
                collapsing.method = collapsing.method, adjust.p.values = adjust.p.values, 
                statistics.type = statistics.type)
        }
        names(res) = nmres
        if (verbose == TRUE) {
            cat("\n")
        }
    }
    if (is.null(model.formula)) {
        names(res) = apply(contrast[ind, , drop = FALSE], 1, 
            function(x) paste(x, collapse = " vs "))
    }
    if (is.null(gsmap)) {
        gsmap = "NULL"
    }
    attr(res, "map") = gsmap
    attr(res, "cov") = covariates
    attr(res, "data") = X
    attr(res, "parameters") = list(X = X, I = I, contrast = contrast, 
        covariates = covariates, model.formula = model.formula, 
        model.type = model.type, ps2gs = ps2gs, gsmap = gsmap, 
        padj.method = padj.method, data.type = data.type, collapsing.method = collapsing.method, 
        include.ordinary.pvalues = include.ordinary.pvalues)
    attr(res, "SessionInfo") = sessionInfo()
    attr(res, "lmfit") = keepFit
    tm = paste0(format(Sys.time(), "%a_%b_%d_%Y_%Hh%Mm%Ss"), 
        "-", paste(sample(c(letters[1:16], LETTERS[1:26]), 25, 
            replace = TRUE), collapse = ""))
    attr(res, "original.names") = paste0(names(res), " (", tm, 
        ")")
    attr(res, "subset") = 1:length(res)
    return(res)
}
