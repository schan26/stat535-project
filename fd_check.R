# Doing a check by comparing our gradient with centered finite differentiation 
# (lecture 1; using centered version improves error a bit more compared to not centered)
# I'm checking both our test functions g1 and g2 with 5 datasets and 5 random bandwidth (just want to see if it matches well for a few sets) for M = 0, 1, 3, 5

source("core.R")
source("grad.R")

set.seed(143)
eps <- 1e-5 # Computed in notes
fun_list <- list(g1 = list(g = g1, sigma = sqrt(0.5)), g2 = list(g = g2, sigma = 0.2))

# Central finite difference
fd_grad <- function(theta, M, X, y){
    sapply(1:3, function(a){
        e <- replace(numeric(3), a, eps)
        (loocv(theta + e, M, X, y) - loocv(theta - e, M, X, y)) / (2 * eps)
    })
}

out <- NULL

# Loop for each function, dataset
for(fname in names(fun_list)){
    for(r in 1:5){
        d <- gen_data(100, fun_list[[fname]]$g, fun_list[[fname]]$sigma)

        for(b in 1:5){

            # random widths (approx h11, h22)
            h <- runif(2, 0.05, 0.3)

            # Gives random theta_2
            theta <- c(log(sqrt(h[1])), runif(1, -0.1, 0.1), log(sqrt(h[2])))

            for(M in c(0, 1, 3, 5)){
                exact <- loocv_grad(theta, M, d$x_data, d$y)
                fd <- fd_grad(theta, M, d$x_data, d$y)

                # Relative error comp
                out <- rbind(out, data.frame(fun = fname, M = M, rel = sqrt(sum((exact - fd)^2))/ sqrt(sum(fd^2))))
            }
        }
    }
}

# worst case over 5 datasets x 5 bandwidths, per function and M
print(aggregate(rel ~ fun + M, out, function(v) signif(max(v), 2)))
cat("Overall worst relative difference:", signif(max(out$rel), 2), "over", nrow(out), "checks\n")