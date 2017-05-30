#' Solving an LRE model according to Blanchard and Kahn
#'
#' @param A Coefficient matrix for the step function
#' @param n Number of predetermined variables
#'
#' @return Two functions, decision rule g() and a motion function h()
#' @export
#'
#' @examples
lre_auto_bk <- function(A, n) {
  sch <- QZ::qz(A)
  sch <- QZ::qz.dtrsen(sch$T, sch$Q, abs(sch$W) <= 1)
  Q1S <- sch$Q[1:n,1:n]
  Q2S <- sch$Q[(n+1):(nrow(sch$Q)),1:n]
  g <- function(x0) Q2S %*% solve(Q1S) %*% x0
  h <- function(x0) A[1:n,1:n] %*% x0 + A[1:n,(n+1):ncol(A)] %*% Q2S %*% solve(Q1S) %*% x0
  return(list(g,h))
}

##

#' Solving an LRE model according to the Klein method with use of QZ decomposition
#'
#' @param A Coefficient matrix on the previous step of the difference equation
#' @param E Coefficient matrix of the current(next) step 
#' @param n Number of predeterminate parameters
#'
#' @return Two functions g() and h() for a decision rule and law of motion respectively
#' @export
#'
#' @examples
lre_auto_klein <- function(A, E, n) {
  
  sch <- QZ::ordqz(E,A, keyword = "udo")
  Z21 <- sch$Z[n:(n+1), 1:n] ## not certain about this
  Z11 <- sch$Z[1:n, 1:n]
  S11 <- sch$S[1:n, 1:n]
  T11 <- sch$T[1:n, 1:n]
  ST <- solve(S11) %*% T11
  g <- function(x0) Z21 %*% solve(Z11) %*% x0
  h <- function(x0) Z11 %*% ST %*% solve(Z11) %*% x0
  return(list(g,h))
  
}

##

#' A combination of the Klein and Blanchard, Kahn methods. 
#'
#' @param A Coefficient matrix on x(t)
#' @param E Coefficient matrix on x(t+1), if not given, the function evaluates the BK method
#' @param n number of predetermined variables
#'
#' @return Motion and decision rule functions according to either the Blanchard Kahn approach or through the Klein method
#' @export
#'
#' @examples
lre_auto <- function(A, E = NULL, n) {
  if (is.null(E)) {
    lre_auto_bk(A, n)
  } else {
    lre_auto_klein(A, E, n)
  }
}
