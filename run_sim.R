# Simulation
# Output: data frame (one row per optimizer run), list (paths) for fig 2, list (tasks) for every data and start point

library(parallel)
source("core.R")
source("grad.R")
source("gd.R")

# Set up
seed <- 143
n <- 100
n_rep <- 10 # Number of datasets to test
n_start <- 10 # Random starts per dataset
M_vals <- c(0, 1, 3, 5)
methods <- c("NM", "GD", "BFGS")
fun_list <- list(g1 = list(g = g1, sigma = sqrt(0.5)), g2 = list(g = g2, sigma = 0.2))
n_cores <- max(1, detectCores() - 1) # set = 1 if not parallel run

# Evaluation grid
grid1 <- seq(0, 1, length.out = 50) # Grid on [0,1]
Xg <- as.matrix(expand.grid(x1 = grid1, x2 = grid1)) # Gives every combo of x1 and x2
margin <- 0.15 # To know if interior/boundary (will update if run gives big h)
is_bd <- Xg[, 1] < margin | Xg[, 1] > 1 - margin | Xg[, 2] < margin | Xg[, 2] > 1 - margin # Checks if it's boundary (TRUE if within 0.15 of edge)

# Estimates error (uses ISE)
est_error <- function(ghat, gtrue){
    c(
        interior = mean((ghat-gtrue)[!is_bd]^2), # interior points are [0.15, 0.85]^2
        boundary = mean((ghat-gtrue)[is_bd]^2),
        overall = mean((ghat-gtrue)^2))
}

# Counter that counts the LOOCV and gradient LOOCV
# Every call from every method goes through f and g so the costs are counted on the same "scale" for all methods
make_counter <- function(M, X, y){
    env <- new.env()
    env$nf <- 0 # Number of LOOCV calls
    env$ng <- 0 # Number of gradient calls
    env$path <- numeric(0) # Stores every LOOCV value in the order it was computed
    env$path_cost <- numeric(0) # Cost 'so far' at each LOOCV evaluation (x-axis for figure 2)

    # LOOCV at each evaluation (LOOCV(theta; M) for a specific dataset)
    f <- function(theta){
        env$nf <- env$nf + 1
        v <- loocv(theta, M, X, y)
        env$path <- c(env$path, v)
        env$path_cost <- c(env$path_cost, env$nf + 6 * env$ng) 
        v
    }

    # Gradient of LOOCV function
    g <- function(theta){
        env$ng <- env$ng + 1
        loocv_grad(theta, M, X, y)
    }

    list(f = f, g = g, env = env)
}

# Per optimizer per starting point per M; th0: starting values for theta
run_one <- function(method, th0, M, X, y){
    cnt <- make_counter(M, X, y) # Creates a new counter per run

    # tryCatch converts singular matrix/issues to NULL so code doesn't stop
    fit <- tryCatch(switch(method,
    NM = optim(th0, cnt$f, method = "Nelder-Mead"), # LOOCV only 
    BFGS = optim(th0, cnt$f, cnt$g, method = "BFGS", control = list(maxit = 500)), # 500 to match other methods' cap/stop criteria
    GD = gd_backtrack(th0, cnt$f, cnt$g, alpha0 = 10)), # alpha0 starting value for decrease
    error = function(e) NULL)
    
    # NULL fit means optimizer crashed/failed or at Inf (bandwidth too tiny)
    if(is.null(fit) || !is.finite(fit$value)){
        return(NULL)
    }
    
    # Want to see if the run stopped because of the optimizer's own stopping criteria and not the iteration max 500
    # loocv is called again here so check isn't counted in count
    # ||gradient of f(theta hat)|| < 10^-6 for gradient descent; BFGS and NM (0 = converged by own rule, 1 = hit iteration cap, 10 = NM only (simplex collapsed))
    converged <- if(method == "GD") sqrt(sum(loocv_grad(fit$par, M, X, y)^2)) < 1e-6 
        else fit$convergence == 0 # Proper stop/convergence
    
    # Cost makes the counting equivalent for all methods since some methods don't use gradient and gradient is expensive
    # 6: finite differentiation with 3 parameters needs about 6 LOOCV calls (2 x 3)
    # cummin = best LOOCV 'so far'
    list(theta = fit$par, value = fit$value, converged = converged, nf = cnt$env$nf, ng = cnt$env$ng, cost = cnt$env$nf + 6 * cnt$env$ng, path = cummin(cnt$env$path), path_cost = cnt$env$path_cost)
}

# Generate every data and their starts
set.seed(seed)
tasks <- list()

# for each test function g1, g2
for (fname in names(fun_list)){

    # for each dataset we are running this on
    for(r in 1:n_rep){
        d <- gen_data(n, fun_list[[fname]]$g, fun_list[[fname]]$sigma)
        starts <- replicate(n_start, {
            # Same starting values for each method and M
            h <- runif(2, 0.05, 0.3)
            theta_from_h(h[1], h[2])
        }, simplify = FALSE)
        
        tasks[[length(tasks) + 1]] <- list(fname = fname, r = r, d = d, starts = starts)
    }
}

# Each data set and one M at a time with each pair using the same data and starting points from above
# start at largest M first (parallel stuffs)
jobs <- list()
for (M in sort(M_vals, decreasing = TRUE)){
    for(k in seq_along(tasks)){
        
        jobs[[length(jobs) + 1]] <- list(task = tasks[[k]], M = M)
    }
}

# All optimizers run for one job (3 methods and n_start starts)
run_job <- function(job){
    task <- job$task # list with function, dataset number, data generated, starts
    M <- job$M
    d <- task$d
    gtrue <- fun_list[[task$fname]]$g(Xg[, 1], Xg[, 2]) # true function on the grid
    rows <- list()
    paths <- list()

    for(meth in methods){
        for(s in seq_len(n_start)){
            out <- run_one(meth, task$starts[[s]], M, d$x_data, d$y) # out is a long list of stuffs that has everything
            row <- data.frame(fun = task$fname, rep = task$r, M = M, method = meth, start = s, failed = is.null(out), converged = NA, loocv = NA, cost = NA, nf = NA, ng = NA, h11 = NA, h22 = NA, h12 = NA, interior = NA, boundary = NA, overall = NA)

            # if not NULL we fill up row that's empty at first above (so record results)
            if(!is.null(out)){
                row$converged <- out$converged
                row[, c("loocv", "cost", "nf", "ng")] <- c(out$value, out$cost, out$nf, out$ng)
                row[, c("h11", "h22", "h12")] <- gen_A(out$theta) # chosen bandwidth from theta hat
                row[ , c("interior", "boundary", "overall")] <- est_error(predict_grid(out$theta, M, d$x_data, d$y, Xg), gtrue) # error of ghat on the grid

                # r = 1 is one dataset and s = 1 is one set of start values
                if(task$r == 1 && s == 1){
                    paths[[paste(task$fname, M, meth)]] <- list(value = out$path, cost = out$path_cost) # Figure 2 uses this stuff
                }
            }
            rows[[length(rows) + 1]] <- row
        }
    }
    list(res = do.call(rbind, rows), paths = paths)
}

# Run everything in parallel
t0 <- Sys.time()

cat("Running", length(jobs), "jobs (", length(tasks), "datasets x", length(M_vals), "values of M), ", length(methods)* n_start, "runs each, on", n_cores, "cores...\n")

# NEED THIS FOR IT TO RUN PARALLEL
if (n_cores > 1) {
  cl <- makeCluster(n_cores)                          # separate R sessions; works on macOS and Windows
  wd <- getwd()
  clusterExport(cl, c("wd", "fun_list", "M_vals", "methods", "n_start", "Xg", "is_bd",
                      "est_error", "make_counter", "run_one"))
  clusterEvalQ(cl, { setwd(wd); source("core.R"); source("grad.R"); source("gd.R"); NULL })
  # parLapplyLB hands out jobs one at a time as cores free up (parLapply would split them up front)
  out <- tryCatch(parLapplyLB(cl, jobs, run_job), finally = stopCluster(cl))
} else {
  out <- lapply(jobs, run_job)
}

# Save results and combine
res <- do.call(rbind, lapply(out, '[[', "res"))
res <- res[order(res$fun, res$rep, res$M, match(res$method, methods), res$start), ]

rownames(res) <- NULL
paths <- do.call(c, lapply(out, '[[', "paths"))

cat("Runs: ", nrow(res), " failed:", sum(res$failed), " time:", round(as.numeric(difftime(Sys.time(),t0, units = "mins")), 1), "min\n")

# Succes if within 1e-4 of the best LOOCV found on that dataset and M
best <- ave(res$loocv, res$fun, res$rep, res$M, FUN = function(v) min(v, na.rm = TRUE))
res$success <- !res$failed & res$loocv <= best * (1 + 1e-4)

save(res, paths, tasks, file = "sim_results.RData")