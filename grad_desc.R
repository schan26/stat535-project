# gradient descent
grad_desc <- function(theta, f, g, alpha0 = 1, rho = 0.5, tau = 1e-4, 
                        tol = 1e-6, maxit=500) {
    fx <- f(theta); num_f <- 1; num_g <- 0; path <- fx

    for (k in 1:maxit) {
        grad <- g(theta); num_g <- num_g + 1
        # gradient ~ 0: stop (reached minimum)
        if (sqrt(sum(grad^2)) < tol) 
            break 

        # steepest line search (negative gradient)
        s <- -grad; 

        # backtracking line search of step size
        alpha <- alpha0                    
        repeat {  # check if need to decay step size                                     
            fnew <- f(theta + alpha * s); num_f <- num_f + 1
            if (fnew <= fx + alpha * tau * sum(s * grad))  # sufficient decrease
                break   
            alpha <- rho * alpha
            if (alpha < 1e-12) # step size decayed to very small
                return(list(par = theta, value = fx, iter = k, nf = num_f, ng = num_g, path = path))
        }
        
        # update theta
        theta <- theta + alpha * s
        fx <- fnew; path <- c(path, fx) # record path
    }
    list(par = theta, value = fx, iter = k, nf = num_f, ng = num_g, path = path)
}
