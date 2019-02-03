
#' Calculating average sojourn time in each state

#' @description Calculating expectation of sojourn times in states (currently two states)
#' for the observed time and for given initial state, using eigenvalues and eigenvectors.
#' @details
#' Calculating expectation of sojourn times in states 1 and 2 for the observed time (tau_observed) and if initial state is given (ii).
#' Matrix Q is so-called Generator matrix: Q = lambda – Lambda, where lambda is matrix with transition rates from state s_i to state s_j,
#' and Lambda is diagonal matrix with a vector (Lambda_1, …, Lambda_m) on the main diagonal, where m is a number of states of external environment (currently m = 2).
#' Eigenvalues and eigenvectors are used in calculations.
#' @param ii
#' number (scalar), currently 1 or 2
#' @param tau_observed
#' number (scalar), observed time
#' @param Q
#' Matrix (m x m, currently m = 2, see Details
#' @return
#' Vector of average sojourn times in each state.
#' Vector components in total should give observation time (tau_observed).
#' @export

#' @examples
#' lambda <- matrix(c(0, 0.33, 0.45, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' m <- nrow(lambda)
#' ld <- as.matrix(rowSums(lambda))
#' Lambda <- diag(as.vector(ld))
#' Generator <- t(lambda) - Lambda
#' Aver_soj_time(1,10,Generator)
#' \donttest{Aver_soj_time(2,5,Generator)}

Aver_soj_time <- function(ii,tau_observed,Q) {
  m <- nrow(Q)
  eigen_vectors <- eigen(Q)$vector
  eigenvals <- as.matrix(eigen(Q)$values)
  b <- as.matrix(solve(eigen_vectors)[,ii])
  R <- array(0,dim = c(m,1))
  for (i in 1:m) {
    if (abs(eigenvals[i,])>1e-10) {
      R <- R + b[i,] * ( (exp(tau_observed*eigenvals[i,])-1) / eigenvals[i,]) * eigen_vectors[,i] }
    else
    { R <- R + b[i,] * eigen_vectors[,i] * tau_observed }
  }
  R
}

#  ---------------Preparing data-------------------------

#' Preparing data for parameter estimation procedure
#' @description
#' Regressors matrix formation taking into account observation times and initial states.
#' Kronecker product is used.
#' @details
#' Function calculates the following expression
#' addImage(doc, "matrix.png")
#' where vector of average sojourn times in each state is calculated using function Aver_soj_time.
#' \figure{matrix.png}{Fig.1}
#' @param tGiven
#' Vector n
#' @param initState
#' Vector n
#' @param X
#' Matrix (n x k+1), k – number of regressors (plus intercept)
#' @param lambda
#' Matrix (m x m), m – number of states
#' @return
#' Matrix (n x 2(k+1))
#' @export
#'
#' @examples
#' X <- cbind(rep_len(1,10),c(2,5,7,3,1,1,2,2,3,6), c(5,9,1,2,3,2,3,5,2,2))
#' tGiven <- matrix (c(6,4.8,1,2.6,6.4,1.7,2.9,4.4,1.5,3.4), nrow = 10, ncol = 1)
#' initState <- matrix (c(2,1,1,2,2,2,1,1,2,1),nrow = 10, ncol = 1)
#' lambda <- matrix(c(0, 0.33, 0.45, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' Xreg(tGiven, initState, X, lambda)

Xreg <- function(tGiven, initState, X, lambda) {
  n <- nrow(X)
  m <- nrow(lambda)
  k <- ncol(X)
  ld <- as.matrix(rowSums(lambda))
  Lambda <- diag(as.vector(ld))
  Q <- t(lambda) - Lambda
  R <- array(0, dim=c(n, m*k))
  for (r in 1:n) {
    ii <- initState[r,]
    taur <- tGiven[r,]
    A <- t(Aver_soj_time(ii,taur,Q))
    B <- t(t(X)[,r])
    R[r,] <- t(kronecker(A, B))
  }
  R }

#-------------Estimantion of unknown model parameters--------------------


#' Estimantion of unknown Markov-modulated linear regression model parameters using GLSM
#' @description
#' This function is used to fit Markov-modulated linear regression models with two states of external environment.
#' This function estimates Markov-modulated linear regression model parameters, using GLSM.
#' Function uses the algorithm based on eigenvalues and eigenvectors decompositions.
#' @details
#' Function calculates the following expression:
#' addImage(doc, "vecB.png")
#' where vector of average sojourn times in each state t_i is calculated using function Aver_soj_time, t_i is an element of tGiven, x_i is a vector of matrix X.

#' @param tGiven
#' Vector n, n – number of observations
#' @param initState
#' Vector n, n – number of observations
#' @param X
#' Matrix (n x k+1), k – number of regressors (plus intercept)
#' @param Y
#' Vector n, n – number of observations
#' @param lambda
#' Matrix (m x m), m – number of states
#' @param W
#' an optional logical variable indicating should vector of weights be used in the fitting process or not.
#' If TRUE, matrix with weights is used (that is, inverse values to tGiven – observed times).
#' @return
#' Vector of estimated parameters
#' @export
#'
#' @examples
#'
#' lambda <- matrix(c(0, 0.33, 0.45, 0), nrow = 2, ncol = 2, byrow = TRUE)
#' X <- cbind(rep_len(1,10),c(2,5,7,3,1,1,2,2,3,6), c(5,9,1,2,3,2,3,5,2,2))
#' tGiven <- matrix (c(6,4.8,1,2.6,6.4,1.7,2.9,4.4,1.5,3.4), nrow = 10, ncol = 1)
#' Y <- matrix(c (5.7, 9.9, 7.8, 14.5, 8.2, 14.5, 32.2, 3.8, 16.5, 7.7),nrow = 10, ncol = 1)
#' initState <- matrix (c(2,1,1,2,2,2,1,1,2,1),nrow = 10, ncol = 1)
#' B_est2(tGiven,initState,X,Y,lambda,W = 1)


B_est2 <- function(tGiven,initState,X,Y,lambda,W = TRUE) {

  k <- ncol(X)		# number of vars
  n <- nrow(X)		# number of obs
  m <- nrow(lambda)

  if (W == TRUE)
    WF <- diag(as.vector(1/tGiven))
  else if (W == FALSE) WF <- diag(as.vector(array(1,dim=c(n,1))))
  else {		stop()
    geterrmessage() }
 # MX <- Xreg(tau=tau, I=I, X=X)
 # MXX <- t(MX) %*% W %*% MX
 # VY <- t(MX) %*% W %*% Y
 b <- solve(t(Xreg(tGiven, initState, X, lambda)) %*% WF %*% Xreg(tGiven, initState, X, lambda)) %*% (t(Xreg(tGiven, initState, X, lambda)) %*% WF %*% Y)

 return(b)}

  # ------------Future Options-----------------------------------------------

  # ----------- Ordinary estimators without environment----------------------

  # D <- diag(as.vector(1/tGiven))

  # bOrdX_Simple <- function (X,Y) {
  #  R <- solve(t(X) %*% X) %*% t(X) %*% Y
  #  R }
  # res1 <- bOrdX_Simple(X,Y)

  # bOrdX <- function(tau, X,Y) {
  #  R <- solve ( t(X) %*% diag(as.vector(1/tau)) %*% X ) %*% t(X) %*% diag(as.vector(1/tGiven)) %*% Y
  #  R }
  # res2 <- bOrdX(tGiven,X,Y)

  # ------------------- RANDOM ENVIRONMENT-------------------------------
  # n <- nrow(Y)
  # m <- nrow(lambda)
  # ld <- as.matrix(rowSums(lambda))
  # Lambda <- diag(as.vector(ld))
  # Generator <- t(lambda) - Lambda

  # ev <- eigen(Generator)
  # eigenvals <- as.vector(ev$values)
  # eigen_vectors <- ev$vectors

  # Generator %*% eigen_vectors[,1]
  # Generator %*% eigen_vectors[,2]

  # eigenvals[1]*eigen_vectors[,1]
  # eigenvals[2]*eigen_vectors[,2]

  #---------------Diagonal matrix---------------------------------
  # Diagonal matrix, exp(eigenval*t) on diagonal

  # Dhi <- function(t, Q) {
  #  Dm <- diag(as.vector(exp(t*as.vector(eigen(Q)$values))))
  #  Dm }

  # Dhi(1, Generator)

  #  --------Calculating transition probability matrix for time delta t (dt)--------
  # ---------Future Option----------------------------------------------------------
  # Transition.pr <- function(t, Q) {
  #  vectors <- eigen(Q)$vectors
  #  pp <- vectors %*% Dhi(t=t, Q=Q) %*% solve(vectors)
  #  pp }

  # ----------Stationary Probabilities--------------------------------------

  # StatPr <- array(0, dim = c(m,1))
  # StatPr[1] <- ld[2]/(ld[1] + ld[2])
  # StatPr[2] <- 1-StatPr[1]
  # StatPr
