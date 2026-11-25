# landscape.R: Fig. 1 (LOOCV contours with optimizer paths) and Table 2 (the landscape in numbers)
# Run AFTER run_sim.R, from the same folder. Uses the datasets saved in sim_results.RData,
# so the landscapes are for exactly the same data as the simulation.
# Output: fig1_landscape.pdf, table2.csv
source("core.R")      # loocv, gen_A, theta_from_h
source("grad.R")      # loocv_grad
source("grad_desc.R")        # gd_backtrack
load("sim_results.RData")   # res, paths, tasks

# ---------- settings ----------
M_show  <- c(0, 5)              # sharpening levels shown in Fig. 1 (columns)
rep_fig <- 1                    # which dataset to draw (same one as Fig. 2)
n_grid  <- 40                   # 40 x 40 grid of bandwidths per panel
h_range <- c(0.02, 0.6)         # range of h11 and h22 shown (log-spaced)
h_start <- c(0.4, 0.4)          # common start for the drawn paths: large widths (the flat side)
methods <- c("NM", "GD", "BFGS")
cols <- c(NM = "#E69F00", GD = "#0072B2", BFGS = "#009E73")   # colour-blind-safe
ltys <- c(NM = 1, GD = 2, BFGS = 4)

# With no tilt (theta2 = 0): h11 = exp(2 theta1), h22 = exp(2 theta3), so theta = log(h) / 2
h_seq  <- exp(seq(log(h_range[1]), log(h_range[2]), length.out = n_grid))
th_seq <- log(h_seq) / 2

get_task <- function(fname, r) tasks[[which(sapply(tasks, function(t) t$fname == fname && t$r == r))]]

# ---------- LOOCV on the (theta1, theta3) grid, tilt fixed at 0 ----------
loocv_grid <- function(M, X, y) {
  z <- matrix(NA, n_grid, n_grid)
  for (i in 1:n_grid) for (j in 1:n_grid) z[i, j] <- loocv(c(th_seq[i], 0, th_seq[j]), M, X, y)
  z[!is.finite(z)] <- NA         # bandwidth too small for the LOOCV formula (P_ii ~ 0)
  z
}

# ---------- best-so-far path of one optimizer run (theta each time LOOCV improves) ----------
trace_run <- function(method, th0, M, X, y) {
  env <- new.env(); env$best <- Inf; env$pts <- matrix(th0, nrow = 1)
  f <- function(t) {
    v <- loocv(t, M, X, y)
    if (is.finite(v) && v < env$best) { env$best <- v; env$pts <- rbind(env$pts, t) }
    v
  }
  g <- function(t) loocv_grad(t, M, X, y)
  tryCatch(switch(method,
    NM   = optim(th0, f, method = "Nelder-Mead"),
    BFGS = optim(th0, f, g, method = "BFGS", control = list(maxit = 500)),
    GD   = grad_desc(th0, f, g)),
    error = function(e) NULL)
  env$pts
}

# ---------- Fig. 1 ----------
th0  <- theta_from_h(h_start[1], h_start[2])
funs <- unique(res$fun)
pdf("fig1_landscape.pdf", width = 9, height = 4.5 * length(funs))
par(mfrow = c(length(funs), length(M_show)), mar = c(4.5, 4.5, 5, 4))
for (fname in funs) {
  d <- get_task(fname, rep_fig)$d
  for (M in M_show) {
    cat("Fig. 1:", fname, "M =", M, "\n")
    z  <- loocv_grid(M, d$x_data, d$y)
    lz <- log10(z)
    lev <- quantile(lz, probs = seq(0.02, 0.98, length.out = 15), na.rm = TRUE)
    contour(th_seq, th_seq, lz, levels = lev, drawlabels = FALSE, col = "grey55",
            xlab = expression(theta[1] == log(h[11]) / 2), ylab = expression(theta[3] == log(h[22]) / 2),
            main = "")
    title(main = paste0(fname, ",  M = ", M), line = 3.2)
    h_ticks <- c(0.02, 0.05, 0.1, 0.2, 0.5)
    axis(3, at = log(h_ticks) / 2, labels = h_ticks, cex.axis = 0.8); mtext(expression(h[11]), 3, 1.9, cex = 0.8)
    axis(4, at = log(h_ticks) / 2, labels = h_ticks, cex.axis = 0.8); mtext(expression(h[22]), 4, 2.2, cex = 0.8)
    k <- which(lz == min(lz, na.rm = TRUE), arr.ind = TRUE)[1, ]
    points(th_seq[k[1]], th_seq[k[2]], pch = 4, cex = 1.6, lwd = 2)            # grid minimum
    for (meth in methods) {
      pts <- trace_run(meth, th0, M, d$x_data, d$y)
      # paths move in 3-D (the tilt theta2 can change); drawn projected onto (theta1, theta3), so an
      # endpoint can differ from the grid minimum, which is the minimum of the tilt = 0 slice only
      lines(pts[, 1], pts[, 3], col = cols[meth], lty = ltys[meth], lwd = 2)
      points(pts[nrow(pts), 1], pts[nrow(pts), 3], col = cols[meth], pch = 19)
    }
    points(th0[1], th0[3], pch = 1, cex = 1.6, lwd = 2)                         # start
    legend("topleft", c(methods, "start", "grid min"), col = c(cols, "black", "black"),
           lty = c(ltys, NA, NA), pch = c(NA, NA, NA, 1, 4), lwd = 2, bty = "n", cex = 0.75, bg = "white")
  }
}
dev.off()

# ---------- Table 2: best bandwidth per M, and how lopsided the valley is ----------
# For every dataset and M: the best theta found by any method/start, then LOOCV with both widths
# halved and doubled (whole H^1/2 scaled by 0.5 or 2), divided by the minimum LOOCV.
rows <- list()
for (t in tasks) {
  for (M in sort(unique(res$M))) {
    sub <- res[res$fun == t$fname & res$rep == t$r & res$M == M & !res$failed, ]
    b <- sub[which.min(sub$loocv), ]
    f_at <- function(c) loocv(theta_from_h(c * b$h11, c * b$h22, c * b$h12), M, t$d$x_data, t$d$y)
    rows[[length(rows) + 1]] <- data.frame(fun = t$fname, rep = t$r, M = M,
      h11 = b$h11, h22 = b$h22, h12 = b$h12, min_loocv = b$loocv,
      ratio_half = f_at(0.5) / b$loocv, ratio_double = f_at(2) / b$loocv)
  }
}
tab2_all <- do.call(rbind, rows)
tab2 <- aggregate(cbind(h11, h22, h12, min_loocv, ratio_half, ratio_double) ~ fun + M,
                  data = tab2_all, FUN = median)            # cells: median over datasets
tab2 <- tab2[order(tab2$fun, tab2$M), ]
num <- sapply(tab2, is.numeric) & names(tab2) != "M"
tab2[num] <- lapply(tab2[num], signif, 3)
print(tab2, row.names = FALSE)
write.csv(tab2, "table2.csv", row.names = FALSE)
cat("Saved fig1_landscape.pdf and table2.csv\n")