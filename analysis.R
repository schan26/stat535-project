# analysis.R: Tables 1 and 3, Fig. 2, appendix tables, and the 0.15-margin check
# Run AFTER run_sim.R, from the same folder. Only reads sim_results.RData (no new LOOCV evaluations).
# Output: table1_<fun>.csv, table3_<fun>.csv, appendix_overall_error_<fun>.csv,
#         appendix_raw_counts.csv, fig2_convergence.pdf
load("sim_results.RData")   # res, paths, tasks

methods <- c("NM", "GD", "BFGS")
M_vals  <- sort(unique(res$M))
funs    <- unique(res$fun)
cols <- c(NM = "#E69F00", GD = "#0072B2", BFGS = "#009E73")   # same colours as Fig. 1
ltys <- c(NM = 1, GD = 2, BFGS = 4)

# ---------- 0. checks before writing ----------
cat("Runs:", nrow(res), "  failed:", sum(res$failed), "\n")
# margin check, per test function: 2 x median selected width vs the 0.15 margin
# (about 0.1-0.2 supports "the margin is about twice the median selected bandwidth")
for (fname in funs) {
  d <- res[res$fun == fname & !res$failed, ]
  cat(fname, ": median h11 =", signif(median(d$h11), 3), "  median h22 =", signif(median(d$h22), 3),
      "  2 x median width =", signif(2 * median(c(d$h11, d$h22)), 3), "\n")
}
cat("\n")

# helper: one row per method, one column per M, from a per-(method, M) function
wide <- function(d, cell) {
  out <- data.frame(method = methods)
  for (M in M_vals) out[[paste0("M=", M)]] <- sapply(methods, function(m) cell(d[d$method == m & d$M == M, ]))
  out
}
pct <- function(x) sprintf("%.0f%%", 100 * mean(x))

# ---------- Table 1: cost, success, converged ----------
# cost: median [IQR] LOOCV-equivalents over completed runs; success and converged: % of ALL runs
# (a failed run counts as neither a success nor converged)
res$conv01 <- !is.na(res$converged) & res$converged
for (fname in funs) {
  d <- res[res$fun == fname, ]
  ok <- d[!d$failed, ]
  cost <- wide(ok, function(x) sprintf("%.0f [%.0f, %.0f]", median(x$cost),
                                       quantile(x$cost, 0.25), quantile(x$cost, 0.75)))
  succ <- wide(d, function(x) pct(x$success))
  conv <- wide(d, function(x) pct(x$conv01))
  cat("Table 1,", fname, "- cost: median [IQR] LOOCV-equivalents\n"); print(cost, row.names = FALSE)
  cat("Table 1,", fname, "- success %\n");   print(succ, row.names = FALSE)
  cat("Table 1,", fname, "- converged %\n"); print(conv, row.names = FALSE); cat("\n")
  names(cost)[-1] <- paste("cost", names(cost)[-1])
  names(succ)[-1] <- paste("success", names(succ)[-1])
  names(conv)[-1] <- paste("converged", names(conv)[-1])
  write.csv(cbind(cost, succ[-1], conv[-1]), paste0("table1_", fname, ".csv"), row.names = FALSE)
}

# ---------- Table 3: interior error, mean (SE) over datasets ----------
# Runs on the same dataset are not independent, so first average over starts within each dataset,
# then take mean and SE over datasets.
err_table <- function(col, fname) {
  ok <- res[res$fun == fname & !res$failed, ]
  per_ds <- aggregate(ok[[col]], by = list(rep = ok$rep, M = ok$M, method = ok$method), FUN = mean)
  wide(per_ds, function(x) sprintf("%.4f (%.4f)", mean(x$x), sd(x$x) / sqrt(nrow(x))))
}
for (fname in funs) {
  t3 <- err_table("interior", fname)
  cat("Table 3,", fname, "- interior error, mean (SE) over datasets\n"); print(t3, row.names = FALSE); cat("\n")
  write.csv(t3, paste0("table3_", fname, ".csv"), row.names = FALSE)
  write.csv(err_table("overall", fname), paste0("appendix_overall_error_", fname, ".csv"), row.names = FALSE)
}

# ---------- Appendix: raw evaluation counts (so readers can use a different cost weight) ----------
ok  <- res[!res$failed, ]
raw <- aggregate(cbind(nf, ng) ~ fun + method + M, data = ok, FUN = median)
raw <- raw[order(raw$fun, raw$method, raw$M), ]
cat("Appendix: median raw counts (nf = LOOCV calls, ng = gradient calls)\n"); print(raw, row.names = FALSE)
write.csv(raw, "appendix_raw_counts.csv", row.names = FALSE)

# ---------- Fig. 2: convergence, log10(best LOOCV so far - f*) vs cost ----------
# One dataset (rep 1, start 1); f* = best LOOCV found by any method/start on that dataset and M.
M_fig <- c(1, 5)
pdf("fig2_convergence.pdf", width = 9, height = 4.5 * length(funs))
par(mfrow = c(length(funs), length(M_fig)), mar = c(4.5, 4.5, 3, 1))
for (fname in funs) for (M in M_fig) {
  fstar <- min(res$loocv[res$fun == fname & res$rep == 1 & res$M == M], na.rm = TRUE)
  curves <- lapply(methods, function(m) paths[[paste(fname, M, m)]])
  names(curves) <- methods
  gap <- function(p) pmax(p$value - fstar, 1e-12)
  xmax <- max(unlist(lapply(curves, function(p) if (is.null(p)) 0 else max(p$cost))))
  yr <- range(log10(unlist(lapply(curves, function(p) if (is.null(p)) NULL else gap(p)))))
  plot(NA, xlim = c(1, xmax), ylim = yr, log = "x", xlab = "cost (LOOCV-equivalents, log scale)",
       ylab = expression(log[10](f - f^"*")), main = paste0(fname, ",  M = ", M))
  for (m in methods) if (!is.null(curves[[m]]))
    lines(curves[[m]]$cost, log10(gap(curves[[m]])), type = "s", col = cols[m], lty = ltys[m], lwd = 2)
  legend("topright", methods, col = cols, lty = ltys, lwd = 2, bty = "n")
}
dev.off()
cat("\nSaved table1_*.csv, table3_*.csv, appendix_*.csv and fig2_convergence.pdf\n")