SimulateMixture <- function(N, w, mu, sigma, nu=4, lambda)
{
    # Number of clusters
    K <- length(w)
    if (K==1) {
        mu <- matrix(mu, 1)
        sigma <- array(sigma, c(1, ncol(mu), ncol(mu)))
    } else if (length(mu)==K) {
        mu <- matrix(mu, K, 1)
        sigma <- array(sigma, c(K, 1, 1))
    }
	# dimension
    py <- ncol(mu)
	y <- matrix(0, N, py)
    nu <- rep(nu,K)
    if (!missing(lambda))
        lambda <- rep(lambda,K)

    label <- sample(1:K, N, replace=T, prob=w)
    count <- table(c(label,1:K))-1
    for (k in 1:K) if (count[k]>0) {    
        y[label==k,] <- rmt(count[k], mu[k,], sigma[k,,], nu[k])
        if (!missing(lambda)) y[label==k,] <- rbox(y[label==k,], lambda[k])    
    }
    y
}
