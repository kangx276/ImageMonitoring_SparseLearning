# Simulate correlated noise.

n <- 16
xgrid <- 1:n
ygrid <- xgrid

RM <- diag(1, nrow = n^2, ncol = n^2)
rho <- 0.8
for (i in 1:n) {
  for (j in 1:n) {
    k1 <- (i - 1) * n + j
    for (i1 in 1:n) {
      for (j1 in 1:n) {
        k2 <- (i1 - 1) * n + j1
        dd <- abs(i - i1)^2 + abs(j - j1)^2
        RM[k1, k2] <- rho^{dd}
        # if ((dd > 0) & (dd <= 10)) {
        #   RM[k1, k2] <- rho^{dd}
        # }
      }
    }
  }
}

sigma <- 0.1 
RMC <- t(chol(RM))
max(abs(RMC %*% t(RMC) - RM))
CovMC <- sigma * RMC


# Wavelet coefficients and their corr matrix

dyn.load("../numerical/imspc.so")
library(parallel)

wav1d_dm <- function(n){
  J <- as.integer(log(n, base = 2))
  out <- .C("wav1d_dm", n = as.integer(n), J = as.integer(J), dm = double(n * 2^J))
  return(matrix(out$dm, nrow = n, byrow = FALSE))
}

set.seed(42)
epsilon <- matrix(rnorm(n^2), nrow = n, byrow = FALSE)
epsvec <- c(t(epsilon)) # vectorize row by row

dm1 <- wav1d_dm(n = n)
dm2 <- wav1d_dm(n = n)
DM <- matrix(0, nrow = nrow(dm1) * nrow(dm2), ncol = ncol(dm1) * ncol(dm2))
for (p1 in 1:ncol(dm1)) {
  for (p2 in 1:ncol(dm2)) {
    p <- (p1 - 1) * ncol(dm2) + p2
    for (i1 in 1:nrow(dm1)) {
      for (i2 in 1:nrow(dm2)) {
        i <- (i1 - 1) * nrow(dm2) + i2
        DM[i, p] <- dm1[i1, p1] * dm2[i2, p2]
      }
    }
  }
}

betahat <- drop(t(DM) %*% epsvec)/(nrow(dm1) * nrow(dm2))

wav2d_coef <- function(Z){
  n1 <- nrow(Z)
  n2 <- ncol(Z)
  J1 <- as.integer(log(n1, base = 2))
  J2 <- as.integer(log(n2, base = 2))
  p <- as.integer(2^{J1 + J2})
  out <- .C("wav2d_coef", n1 = n1, n2 = n2, Z = as.double(t(Z)), J1 = J1, J2 = J2, wcoef = double(p))
  return(list(J1 = out$J1, J2 = out$J2, wcoef = out$wcoef))
}

betahat / wav2d_coef(Z = epsilon)$wcoef # checking that it is the same as the 2-D implementation.

## Correlated noise

n == nrow(dm1)

theoryVal <- t(DM) %*% (sigma^2 * RM) %*% DM / (n^4)
theoryValR <- cov2cor(theoryVal)

## Sample cov matrix

set.seed(4200)
N <- 10^3
EpsMat <- matrix(rnorm(n^2 * N), nrow = n^2, byrow = FALSE)
ErrMat <- CovMC %*% EpsMat
BetaHatMat <- t(DM) %*% ErrMat / (n^2)
EmprVal <- cov(t(BetaHatMat))

max(abs(EmprVal - theoryVal))

## Threshold the sample cov matrix. J. Fan (2015), Econometric Journal

dgvec <- diag(EmprVal)
EmprValR <- cov2cor(EmprVal)
max(abs(EmprValR - theoryValR))
sum((theoryValR - EmprValR)^2)

C <- 1
thresh <- C * sqrt(log(nrow(EmprValR)) / N)
EmprValR.thresh <- sign(EmprValR) * 0.5 * (abs(EmprValR) - thresh +  abs(abs(EmprValR) - thresh)) # soft thresholding
diag(EmprValR.thresh) <- 1
max(abs(EmprValR.thresh - theoryValR))

eivals <- eigen(EmprValR.thresh, symmetric = TRUE, only.values = TRUE)$values

## Search the minimum threshold that ensures positive definiteness

EVmin <- function(CorrMat, nobs, Cseq = seq(0, 3, by = 0.1)) {
  thresh.seq <- Cseq * sqrt(log(nrow(CorrMat)) / nobs)
  EVmin.seq <- rep(NA, length(thresh.seq))
  for (k in seq_along(thresh.seq)) {
    thresh <- thresh.seq[k]
    CorrMat.thresh <- sign(CorrMat) * 0.5 * (abs(CorrMat) - thresh +  abs(abs(CorrMat) - thresh)) # soft thresholding
    diag(CorrMat.thresh) <- 1
    EVmin.seq[k] <- min(eigen(CorrMat.thresh, symmetric = TRUE, only.values = TRUE)$values)
  }
  return(EVmin.seq)
}

## Search the optimal threshold that ensures positive definiteness

Frob <- function(CorrMat, nobs, Cpos, trueMat) {
  thresh.seq <- Cpos * sqrt(log(nrow(CorrMat)) / nobs)
  Frob.seq <- rep(NA, length(thresh.seq))
  for (k in seq_along(thresh.seq)) {
    thresh <- thresh.seq[k]
    CorrMat.thresh <- sign(CorrMat) * 0.5 * (abs(CorrMat) - thresh +  abs(abs(CorrMat) - thresh)) # soft thresholding
    diag(CorrMat.thresh) <- 1
    Frob.seq[k] = sum((CorrMat.thresh - trueMat)^2)
  }
  return(Frob.seq)
}

EVmin(CorrMat = EmprValR, nobs = N, Cseq = seq(3, 5, by = 0.1))
Cpos <- seq(3.6, 5, by = 0.1)
Frob(CorrMat = EmprValR, nobs = N, Cpos = Cpos, trueMat = theoryValR)

## Compute AL estimates

### Generate new observation
set.seed(42)
epsilon <- matrix(rnorm(n^2), nrow = n, byrow = FALSE)
epsvec <- c(t(epsilon)) # vectorize row by row
errvec <- CovMC %*% epsvec
betahat <- drop(t(DM) %*% errvec)/(nrow(dm1) * nrow(dm2)) 

### Turn AL setup into Lasso setup
L <- t(chol(theoryVal))
Lam <- diag(abs(betahat))
eta <- forwardsolve(l = L, x = betahat)
Xi <- forwardsolve(l = L, x = Lam)

### LARS algorithm
library(lars)

lastIndx <- function(k, mvec) {
  return(max(which(mvec == k)))
}

m.last <- function(coef.path, q, eps = 10^{-12}) {
  n.act <- apply(coef.path, 1, FUN = function(x) sum(abs(x) > eps))
  stopifnot(sum(abs(sort(unique(n.act))[1:(q + 1)] - c(0:q))) < eps)
  indx.last <- mapply(FUN = lastIndx, k = 0:q, MoreArgs = list(mvec = n.act))
  return(coef.path[indx.last, ])
}

m.last <- function(coef.path, eps = 10^{-12}) {
  n.act <- apply(coef.path, 1, FUN = function(x) sum(abs(x) > eps))
  indx.last <- mapply(FUN = lastIndx, k = sort(unique(n.act)), MoreArgs = list(mvec = n.act))
  return(coef.path[indx.last, ])
}

W1.al <- function(betaOLS, L.chol, q, eps = .Machine$double.eps) {
  stopifnot(min(abs(betaOLS)) >= eps)
  Lam <- diag(abs(betaOLS))
  eta <- forwardsolve(l = L.chol, x = betaOLS)
  Xi <- forwardsolve(l = L.chol, x = Lam)
  lars.path <- coefficients(lars(x = Xi, y = eta, type = "lasso", normalize = FALSE, intercept = FALSE))
  #stopifnot(nrow(lars.path) >= (q + 1))
  #AL.path <- (m.last(coef.path = lars.path, q = q) %*% Lam)[2:(q + 1), ]
  n.act <- apply(lars.path, 1, FUN = function(x) sum(abs(x) > eps))
  q.real <- sort(unique(n.act[n.act <= q]))[-1] # realized q's in 1:q
  indx.last <- mapply(FUN = lastIndx, k = q.real, MoreArgs = list(mvec = n.act))
  AL.path <- lars.path[indx.last, ] %*% Lam
  xi.path <- forwardsolve(l = L.chol, x = t(AL.path))
  num <- (drop(t(eta) %*% xi.path))^2
  den <- colSums((xi.path)^2)
  out <- rep(NA, q)
  stopifnot(length(num) == length(den))
  stopifnot(length(num) == length(q.real))
  out[q.real] <- num/den
  return(out)
}

head(W1.al(betaOLS = betahat, L.chol = L, q = ncol(DM) - 1))

## SCAD: coordinate descent algorithm

L.inv <- forwardsolve(l = L, x = diag(1, nrow = nrow(L)))
max(abs(L.inv %*% L - diag(1, nrow = n^2, ncol = n^2)))

library(ncvreg)

scadfit <- function(X, y, lambda, eps = 1e-05) {
  return(ncvfit(X = X, y = y, penalty = "SCAD", lambda = lambda, eps = eps, max.iter = 10000)$beta)
}

cda <- function(X, y, lamseq = seq(10, 0.01, by = 0.01), eps = 1e-05) {
  coef.path <- mapply(FUN = scadfit, lambda = lamseq, MoreArgs = list(X = X, y = y, eps = eps))
  return(coef.path)
}

### set up lambda sequence
lambda.max <- max(crossprod(L.inv, eta)) / nrow(L.inv)
#lambda.max <- 2e06
nlambda <- 100
lambda_vec <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
lambda_vec <- seq(3, 1, by = -0.001)

nsims <- 100
succ <- rep(NA, nsims)
for (k in seq_len(nsims)) {
  epsilon <- matrix(rnorm(n^2), nrow = n, byrow = FALSE)
  epsvec <- c(t(epsilon)) # vectorize row by row
  errvec <- CovMC %*% epsvec
  betahat <- drop(t(DM) %*% errvec)/(nrow(dm1) * nrow(dm2)) 
  eta <- forwardsolve(l = L, x = betahat)
  qs <- sort(unique(apply(t(cda(X = L.inv, y = eta, lamseq = lambda_vec)), 1, FUN = function(x) sum(abs(x) > .Machine$double.eps))))
  #print(qs)
  succ[k] <- all(abs(qs[qs > 0][1:10] - 1:10) <= .Machine$double.eps)
}
sum(succ)

W1.scad <- function(betaOLS, L.chol, q, lamseq, eps = .Machine$double.eps) {
  #stopifnot(min(abs(betaOLS)) >= eps)
  #Lam <- diag(abs(betaOLS))
  eta <- forwardsolve(l = L.chol, x = betaOLS)
  Xi <- forwardsolve(l = L.chol, x = diag(1, nrow = length(betaOLS)))
  scad.path <- t(cda(X = Xi, y = eta, lamseq = lamseq))
  n.act <- apply(scad.path, 1, FUN = function(x) sum(abs(x) > eps))
  q.real <- sort(unique(n.act[n.act <= q]))
  q.real <- q.real[q.real > 0]# realized q's in 1:q
  indx.last <- mapply(FUN = lastIndx, k = q.real, MoreArgs = list(mvec = n.act))
  SCAD.path <- scad.path[indx.last, ]# %*% Lam
  xi.path <- forwardsolve(l = L.chol, x = t(SCAD.path))
  num <- (drop(t(eta) %*% xi.path))^2
  den <- colSums((xi.path)^2)
  out <- rep(NA, q)
  stopifnot(length(num) == length(den))
  stopifnot(length(num) == length(q.real))
  out[q.real] <- num/den
  return(out)
}

scad1 <- W1.scad(betaOLS = betahat, L.chol = L, q = ncol(DM) - 1, lamseq = lambda_vec)
length(scad1) - sum(is.na(scad1))
sum(qs <= 255)
head(scad1)

## IC statistics of W1.al

nrep <- 25000
set.seed(314)
ErrMat <- CovMC %*% matrix(rnorm(nrep * n^2), ncol = nrep)
BetaHatMat <- as.data.frame((t(DM) %*% ErrMat) / (nrow(dm1) * nrow(dm2)))
W1.ALmat <- mcmapply(FUN = W1.al, betaOLS = BetaHatMat, MoreArgs = list(L.chol = L, q = ncol(DM) - 1), mc.cores = detectCores())

any(is.na(W1.ALmat))

alo.ic.bar <- rowMeans(W1.ALmat, na.rm = TRUE)
alo.ic.sd <- sqrt(rowMeans(W1.ALmat^2, na.rm = TRUE) - rowMeans(W1.ALmat, na.rm = TRUE)^2)

## Calculate the EWMA of phase II observations

arl0 <- 20#0
n.phase2 <- 10 * arl0
set.seed(442000)
ErrMat.phase2 <- CovMC %*% matrix(rnorm(n.phase2 * n^2), ncol = n.phase2, byrow = FALSE)
BetaHat.phase2 <- (t(DM) %*% ErrMat.phase2) / (nrow(dm1) * nrow(dm2))

la <- 0.2

EWMA <- function(obs.mat, la){
  out <- t(as.matrix(filter(la * t(obs.mat), 1 - la, method = "recursive"))) # convention: the column index is the time index
  return(out)
}

U.phase2 <- EWMA(obs.mat = BetaHat.phase2, la = la)

## Run LARs on the EWMA sequence

W2.al <- function(U.mat, la, L.chol, q) {
  stopifnot(la <= 1 && la > 0)
  t.seq <- seq_len(ncol(U.mat))
  ewma.seq <- (2 - la)/(la * (1 - (1 - la)^{2 * t.seq}))
  W.mat <- mapply(W1.al, betaOLS = as.data.frame(U.mat), MoreArgs = list(L.chol = L.chol, q = q), SIMPLIFY = TRUE)
  out <- t(t(W.mat) * ewma.seq) ## a fast way to multiply each column of W.mat by elements of the EWMA sequence
  return(out)
}

W2.ALmat <- W2.al(U.mat = U.phase2, la = la, L.chol = L, q = 5)
W2.ALmat[, 1:3]

## AL charting statistics: standardize and maximize

Q.al <- function(W.mat, IC.bar, IC.sd) {
  stopifnot(nrow(W.mat) == length(IC.bar))
  stopifnot(length(IC.bar) == length(IC.sd))
  W.centered <- sweep(x = W.mat, MARGIN = 1, STATS = IC.bar, FUN = "-")
  W.sdz <- sweep(x = W.centered, MARGIN = 1, STATS = IC.sd, FUN = "/")
  Q.out <- do.call(pmax.int, c(as.data.frame(t(W.sdz)), na.rm = TRUE))
  return(Q.out)
}

CC.alo <- Q.al(W.mat = W2.ALmat, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd)

## Test the AL chart 

set.seed(442000)
n.ic <- 15
n.oc <- 25
mu.oc <- c(rep(10, 5), rep(0, n^2 - 5))
ErrMat.phase2 <- CovMC %*% cbind(matrix(rnorm(n.ic * n^2), ncol = n.ic, byrow = FALSE), matrix(rnorm(n.oc * n^2), ncol = n.oc, byrow = FALSE) + mu.oc)
BetaHat.phase2 <- (t(DM) %*% ErrMat.phase2) / (nrow(dm1) * nrow(dm2))

U.phase2 <- EWMA(obs.mat = BetaHat.phase2, la = la)
W2.ALmat <- W2.al(U.mat = U.phase2, la = la, L.chol = L, q = 5)
CC.alo <- Q.al(W.mat = W2.ALmat, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd)
plot(CC.alo)

## Average run length

SimQ.al <- function(nobs, ErrRM.chol, wav2DM, L.chol, la, q, IC.bar, IC.sd) { 
  stopifnot(nrow(ErrRM.chol) == ncol(ErrRM.chol))
  stopifnot(nrow(L.chol) == ncol(L.chol))
  stopifnot(nrow(L.chol) == ncol(wav2DM))
  stopifnot(nrow(ErrRM.chol) == nrow(wav2DM))
  stopifnot(la > 0 && la <= 1)
  stopifnot(q == length(IC.bar) && q == length(IC.sd))
  N <- nrow(ErrRM.chol)
  ErrMat <- ErrRM.chol %*% matrix(rnorm(nobs * N), ncol = nobs, byrow = FALSE)
  BetaHat <- t(wav2DM) %*% ErrMat / N
  Umat <- EWMA(obs.mat = BetaHat, la = la)
  W2mat <- W2.al(U.mat = Umat, la = la, L.chol = L.chol, q = q)
  CC <- Q.al(W.mat = W2mat, IC.bar = IC.bar, IC.sd = IC.sd)
  return(CC)
}
simCC <- SimQ.al(nobs = n.phase2, ErrRM.chol = CovMC, wav2DM = DM, L.chol = L, la = la, q = 5, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd)
simCC <- mapply(FUN = SimQ.al, nobs = rep(10 * arl0, 2), 
                MoreArgs = list(ErrRM.chol = CovMC, wav2DM = DM, L.chol = L, la = la, q = 5, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd))

arl0 <- 200
nruns <- 25000
RNGkind("L'Ecuyer-CMRG")
set.seed(2000)
system.time(
  QM.alo <- mcmapply(FUN = SimQ.al, nobs = rep(10 * arl0, nruns), 
                     MoreArgs = list(ErrRM.chol = RMC, wav2DM = DM, L.chol = L, la = la, q = 5, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd),
                     mc.cores = detectCores())
)
invisible(rnorm(1))
dim(QM.alo)
write.csv(QM.alo, file = "QMalo.csv")

## Phase II simulation: AL

f <- function(r) ifelse( r < 0.25, 1, 0)
x <- (0:(n - 1))/n
y <- x;
xy <- as.matrix(expand.grid(x, y))
rr <- sqrt((xy[, 1] - 0.5)^2 + (xy[, 2] - 0.5)^2)
ff <- matrix(f(rr), ncol = n)

image(ff, col = gray(c(0:255)/255))

beta_true <- drop(t(DM) %*% c(t(ff)))/(nrow(dm1) * nrow(dm2))
p <- length(beta_true)
p == n^2
set.seed(314)
Z <- ff + matrix(CovMC %*% rnorm(length(x) * length(y)), ncol = n, byrow = TRUE) # The vectorization of the noise matrix was done row by row. 
image(Z, col = gray(c(0:255)/255))

dyn.load("../numerical/repsim.so")

fault_cross <- function(n1, n2, size) {
  out <- .C("fault_cross", n1 = as.integer(n1), n2 = as.integer(n2), faultsize = as.double(size), dff = double(n1 * n2))$dff
  out <- matrix(out, ncol = n2, byrow = TRUE)
  return(out)
}

fault_square <- function(n1, n2, size) {
  out <- .C("fault_square", n1 = as.integer(n1), n2 = as.integer(n2), faultsize = as.double(size), dff = double(n1 * n2))$dff
  out <- matrix(out, ncol = n2, byrow = TRUE)
  return(out)
}

n_imgs <- 100
tau <- 80
ims_betahat <- matrix(NA, nrow = n^2, ncol = n_imgs)
fs <- 0.4
#dff <- fault_square(n1 = n, n2 = n, size = fs)
dff <- fault_cross(n1 = n, n2 = n, size = fs)
#set.seed(314)
set.seed(42)
for (t in seq_len(n_imgs)) {
  if (t <= tau) {
    ims_betahat[, t] <- drop(t(DM) %*% c(t(ff + matrix(CovMC %*% rnorm(length(x) * length(y)), ncol = n, byrow = TRUE)))) / (nrow(dm1) * nrow(dm2))
  } else {
    ims_betahat[, t] <- drop(t(DM) %*% c(t(ff + dff + matrix(CovMC %*% rnorm(length(x) * length(y)), ncol = n, byrow = TRUE)))) / (nrow(dm1) * nrow(dm2))
  }
}

ims_U <- EWMA(obs.mat = ims_betahat - beta_true, la = la)
ims_W2AL <- W2.al(U.mat = ims_U, la = la, L.chol = L, q = ncol(DM) - 1)
ims_CCAL <- Q.al(W.mat = ims_W2AL, IC.bar = alo.ic.bar, IC.sd = alo.ic.sd)
par(mfrow = c(1, 1))
plot(ims_CCAL[1:n_imgs], type = "b")
lines(x = 0:100, y = rep(4.644, 101), lty = "dashed", col = "red", lwd = 2)
min(which(ims_CCAL >= 4.644))

EI <- eigen(theoryVal, symmetric = TRUE)
EI.val <- EI$values
V <- EI$vectors
max(abs(V %*% diag(EI.val) %*% t(V) - theoryVal))

SGM.irt <- V %*% diag(sqrt(1/EI.val)) %*% t(V)

SGM.inv <- solve(theoryVal)
SGM.rt <- V %*% diag(sqrt(EI.val)) %*% t(V)


W.phase1 <- function(X, S.inv, S.irt, q) {
  Lam <- diag(abs(X))
  x <- S.irt %*% Lam
  y <- S.irt %*% X
  L.path <- coefficients(lars(x = x, y = y, type = "lasso", normalize = FALSE, intercept = FALSE))
  AL.path <- (m.last(L.path) %*% Lam)[2:(q + 1), ]
  num <- c(AL.path %*% S.inv %*% X)^2
  den <- colSums((S.irt %*% t(AL.path))^2)
  return(num/den)
}
W.phase1(X = betahat, S.inv = SGM.inv, S.irt = SGM.irt, q = 5)

W.mat <- mcmapply(FUN = W.phase1, X = BetaHatMat, MoreArgs = list(S.inv = SGM.inv, S.irt = SGM.irt, q = 5), mc.cores = detectCores())

rowMeans(W.mat)
sqrt(rowMeans(W.mat^2) - rowMeans(W.mat)^2)

W.phase2 <- function(U, t, la, S.inv, S.irt, q) {
  stopifnot(la <= 1 && la > 0)
  stopifnot(is.integer(t))
  return((2 - la)/(la * (1 - (1 - la)^{2*t})) * W.phase1(X = U, S.inv = S.inv, S.irt = S.irt, q = q))
  #return(((2 - la)/la) * W.phase1(X = U, S.inv = S.inv, S.irt = S.irt, q = q))
}

W.phase2(U = U.phase2[, 3], t = as.integer(3), la = 0.2, S.inv = SGM.inv, S.irt = SGM.irt, q = 5)

