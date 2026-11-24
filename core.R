# Contains all the functions that will be used

# Create K matrix (Kij)
kernel_matrix <- function(x_eval, x_data, A){
    #Kij: Gaussian kernel of point j (x_data) when estimating evaluation point at i (x_eval)
    # A: vector with h11, h22, h12

    #Build H^1/2 matrix
    H_half <- matrix(c(A[1],A[3],A[3],A[2]),nrow=2)

    # Gets H^-1 matrix
    H_half_inv <- solve(H_half)
    H_inv <- H_half_inv %*% H_half_inv

    # Let uij = xi - xj; build table of values where each row is the evaluation point
    # Sign flip is fine since we square it

    # Coordinate 1 (x1i)
    u1 <- outer(x_eval[, 1], x_data[, 1], FUN = "-")

    # Coordinate 2 (x2i)
    u2 <- outer(x_eval[, 2], x_data[, 2], FUN = "-")

    # Product/Norm: a11u1^2 + 2(a12u1u2) + a22u2^2
    quad <- H_inv[1, 1] * u1^2 + 2 * H_inv[1, 2] * u1 * u2 + H_inv[2, 2] * u2^2

    # Final equation
    exp(-0.5 * quad)
}

# Create Smoother matrix S
S_lc <- function(K){
    K / rowSums(K) # Sij = Kij/di where di = sum over k (K_ik)
}

# Create bandwidth using log-cholesky and optimized theta
gen_A <- function(theta){
    # Log-cholesky: creates a valid bandwidth (positive definite) using any real number

    L <- matrix(c(exp(theta[1]), theta[2], 0, exp(theta[3])), nrow = 2)
    H_half <- L %*% t(L) # H^1/2 = LL^t

    c(H_half[1, 1], H_half[2, 2], H_half[1, 2]) # (h11, h22, h12)
}

# Reverse of gen_A; used for starting points (computes theta given h)
theta_from_h <- function(h11, h22, h12 = 0){
    l11 <- sqrt(h11)
    l21 <- h12 / l11
    l22 <- sqrt(h22 - l21^2)
    c(log(l11), l21, log(l22))
}

# Iterated Data Sharpening function
lc_ids <- function(S, y, M){
    # n x (M + 1) matrix: column M + 1 is the fit at the data points after M iterations

    g_hat_sharp <- matrix(NA, nrow = length(y), ncol = M + 1)

    # At iteration M = 0 (regular lc estimation), (I - S)^0 * y = y
    Py <- y 

    for(m in 0:M){
        Py <- Py - S %*% Py  # (I - S)Py  <=> (I - S)^(m+1) y = y - fit_m
        g_hat_sharp[, m+1] <- y - Py  # fit_m = [I - (I - S)^(m+1)] y = F_m y = y_hat
    }

    g_hat_sharp
}

# Sharpened data
sharpen_y <- function(S, y, M){

    ys <- y # y^(0) = y (starting value of recursions)

    for(m in seq_len(M)){
        ys <- y + (ys - S %*% ys) # y^(m+1) = y^(0) + (I - S)y^(m)
    }

    drop(ys)
}

# Grid of estimates
predict_grid <- function(theta, M, x_data, y, x_grid){

    # IDS estimate at new points x_grid (lc applied to sharpened y^(m))

    A <- gen_A(theta)
    S <- S_lc(kernel_matrix(x_data, x_data, A))
    ys <- sharpen_y(S, y, M)
    Kg <- kernel_matrix(x_grid, x_data, A) # rows = grid points, col = data points
    drop(Kg %*% ys) / rowSums (Kg) # g_hat at iteration M
}

# Objective function LOOCV
loocv <- function(theta, M, x_data, y){
    S <- S_lc(kernel_matrix(x_data, x_data, gen_A(theta)))
    n <- length(y)
    P <- diag(n) # P = (I-S)^(M+1)

    for(m in 0:M){
        P <- P - S %*% P
    } 

    F_M <- diag(n) - P # FM = I-P
    fit <- F_M %*% y # yhat = FM*y
    denom <- diag(P) # Pii

    if(any(!is.finite(fit)) || any(denom <= 1e-8)){
        return(Inf) # if bandwidth is too tiny
    }
    mean(((y - fit) / denom)^2)
}

# Simulation helper function
g1 <- function(x1, x2){
    5 * (x1 - 0.4)^2 + 3 * x1 * x2^2
}

g2 <- function(x1, x2){
    sin(pi * x1) * sin(2 * pi * x2)
}

gen_data <- function(n, g, sigma){
    x_data <- cbind(runif(n), runif(n))
    list(x_data = x_data, y =g(x_data[, 1], x_data[, 2]) + rnorm(n, 0, sigma))
}