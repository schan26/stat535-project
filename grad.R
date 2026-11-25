# Computation of the gradient of LOOCV


## LOOCV = mean(((Py)_i / P_ii)^2), P = (I - S)^{M+1} = R^{M+1}
loocv_grad <- function(theta, M, x_data, y) {
    # number of data points
    n <- length(y)

    # log-cholesky decomposition of H^{1/2}
    L <- matrix(c(exp(theta[1]), theta[2], 0, exp(theta[3])), 2)
    
    # H inverse
    H_half <- L %*% t(L)
    H_half_inv <- solve(H_half)
    H_inv <- H_half_inv %*% t(H_half_inv)

    # K computation
    # Coordinate 1 (x1i)
    u1 <- outer(x_data[, 1], x_data[, 1], FUN = "-")
    # Coordinate 2 (x2i)
    u2 <- outer(x_data[, 2], x_data[, 2], FUN = "-")
    # Product/Norm: a11u1^2 + 2(a12u1u2) + a22u2^2
    quad <- function(A) {
        A[1, 1] * u1^2 + 2 * A[1, 2] * u1 * u2 + A[2, 2] * u2^2
    }
    # K
    K <- exp(-0.5 * quad(H_inv))

    # S = K / d, d_i = rowsum of K_i
    d <- rowSums(K)
    S <- K / d

    # R = I - S
    R <- diag(n) - S

    # P = R^{M+1} computation
    Rp <- vector("list", M + 2)
    Rp[[1]] <- diag(n)                    
    for (k in 1:(M + 1)) Rp[[k + 1]] <- Rp[[k]] %*% R # Rp[[k+1]] = R^k
    P  <- Rp[[M + 2]]

    # r = Py
    r <- drop(P %*% y)
    p <- diag(P) # p_i = P_ii
    

    # error, cannot divide by 0 (bandwidth too small)
    if (any(p <= 1e-8)) return(rep(NA, 3))

    e <- r / p # LOO residuals

    ## compute LOOCV gradient
    loocv_g <- numeric(3)

    # dL / dtheta
    dL <- list(matrix(c(exp(theta[1]),0,0,0),2), matrix(c(0,1,0,0),2),
             matrix(c(0,0,0,exp(theta[3])),2))
    
    # dLOOCV / dtheta_i
    for (i in 1:3) {
        # partial of H^{-1}
        dH_half <- L %*% t(dL[[i]]) + dL[[i]] %*% t(L)
        dH_half_inv <- - H_half_inv %*% dH_half %*% H_half_inv
        dH_inv <- H_half_inv %*% dH_half_inv + dH_half_inv %*% H_half_inv

        # partial of K
        dK <- - K * quad(dH_inv) / 2  #chain rule multiplication

        # partial of S
        dd <- rowSums(dK)
        dS <- (dK - S * dd) / d

        # partial of R
        dR <- -dS

        # partial of P
        dP <- matrix(0, n, n) # initialize dP
        for (k in 0:M) { # product rule
            dP <- dP + Rp[[k + 1]] %*% dR %*% Rp[[M - k + 1]]  
        }

        # partial of p
        dp <- diag(dP)

        # partial of r
        dr <- drop(dP %*% y)

        # partial of g
        loocv_g[i] <- mean(2 * e * (dr - e * dp) / p)
    }

    loocv_g
}
