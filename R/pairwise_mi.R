#' Compute pairwise association matrix wiht mutual information.
#'
#' @export
#'
pairwise_mi <- function(theta,pi) {
    compute_pi = function(i) {
      pi_i = lapply(1:R,function(r) {
        theta[r] * pi[r,i,]
      }) %>% Reduce('+', .)
      pi_i
    }
    compute_pi2 = function(i,j) {
      pi_ij = lapply(1:R,function(r) {
        theta[r] * pi[r,i,] %o% pi[r,j,]
      }) %>% Reduce('+', .)
      pi_ij
    }
    I = function(i,j) {
      A = compute_pi2(i,j)
      B = compute_pi(i) %o% compute_pi(j)
      sum(A * log(A/B))
    }
    H = function(i) {
      A = compute_pi(i)
      -sum(A*log(A))
    }
    M = matrix(nrow = J,ncol = J)
    for(i in 1:J) {
      for(j in 1:J)
        M[i,j] = I(i,j)/sqrt(H(i)*H(j))
    }
    M
}
