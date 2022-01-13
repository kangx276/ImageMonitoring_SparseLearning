rm(list = ls())

p <- 5
s2 <- 0.5^2
SO <- diag(s2, nrow = p)

SO.irt <- diag(1/sqrt(s2), nrow = p)
SO.inv <- diag(1/s2, nrow = p)
SO.rt <- diag(sqrt(s2), nrow = p)

# solve adaptive lasso-type penalized likelihood 

m.last <- function(coef.path, eps = 10^{-12}) {
  n.act <- apply(coef.path, 1, FUN = function(x) sum(abs(x) > eps))
  indx.last <- unique(rank(n.act, ties.method = "max"))
  return(coef.path[indx.last, ])
}

library(lars)

W.phase1 <- function(X, S.inv, S.irt, q) {
  Lam <- diag(abs(X))
  x <- S.irt %*% Lam
  y <- S.irt %*% X
  L.path <- coefficients(lars(x = x, y = y, type = "lasso", normalize = FALSE, intercept = FALSE))
  AL.path <- m.last(L.path %*% Lam)[2:(q + 1), ]
  num <- c(AL.path %*% S.inv %*% X)^2
  den <- colSums((S.irt %*% t(AL.path))^2)
  return(num/den)
}

set.seed(314)
W.phase1(X = c(SO.rt %*% rnorm(p)), S.inv = SO.inv, S.irt = SO.irt, q = 3)

library(parallel)

nrep <- 25000
set.seed(314)
X.mat <- as.data.frame(SO.rt %*% matrix(rnorm(nrep * p), ncol = nrep))
system.time(
  W.mat <- mcmapply(FUN = W.phase1, X = X.mat, MoreArgs = list(S.inv = SO.inv, S.irt = SO.irt, q = 3), mc.cores = detectCores())
)
ic.bar <- rowMeans(W.mat)
ic.sd <- sqrt(rowMeans(W.mat^2) - rowMeans(W.mat)^2)

dyn.load("imspc.so")

ALO_ICstat <- function(n_IC, p, sig2) {
  out <- .C("ALO_ICstat", n_IC = as.integer(n_IC), p = as.integer(p), sig2 = as.double(sig2), IC_bar = double(p), 
            IC_std = double(p))
  ic.stat <- cbind(out$IC_bar, out$IC_std)
  colnames(ic.stat) <- c("mean", "std")
  return(ic.stat)
}

set.seed(314)
ICstat <- ALO_ICstat(n_IC = nrep, p = p, sig2 = s2)

SCAD_phase1 <- function(X, sig2, a = 3.7){
  p = length(X)
  out <- .C("SCAD_phase1", p = p, sig2 = as.double(sig2), a = as.double(a), X = as.double(X), V = double(p))
  return(out$V)
}

pp <- 5
set.seed(34)
ss <- 0.1
xx <- rnorm(pp, sd = ss)

SCAD_phase1(X = xx, sig2 = ss^2)

scad_mu <- function(x, zeta, a = 3.7) {
  #x <- xx[1]
  #zeta <- 0.02
  #a <- 3.7
  sign(x) * ifelse(abs(x) > zeta, abs(x) - zeta, 0) * ifelse(abs(x) <= 2 * zeta, 1, 0) + 
    ((a - 1) * x - sign(x) * a * zeta) /(a - 2) * ifelse((abs(x) > 2 * zeta) & (abs(x) <= a * zeta), 1, 0) + 
    x * ifelse(abs(x) > a * zeta, 1, 0)
}

scad_mu(x = xx, zeta = 0.02)

scad_chart <- function(x, zeta, sig2, q, a = 3.7) {
  p <- length(x)
  if (q < p) {
    shrink <- sort(abs(x))[p - q]
  } else {
    shrink <- 0
  }
  mu <- scad_mu(x = x, zeta = shrink, a = a)
  return(sum(x * mu)^2/(sum(mu^2))/sig2)
}

scad_chart(x = xx, zeta = 0.02, sig2 = ss^2, q = 4, a = 3.7)

SCAD_ICstat <- function(n_IC, p, sig2, a = 3.7) {
  out <- .C("SCAD_ICstat", n_IC = as.integer(n_IC), p = as.integer(p), sig2 = as.double(sig2), a = as.double(a), IC_bar = double(p), 
            IC_std = double(p))
  ic.stat <- cbind(out$IC_bar, out$IC_std)
  colnames(ic.stat) <- c("mean", "std")
  return(ic.stat)
}

set.seed(314)
scadICstat <- SCAD_ICstat(n_IC = nrep, p = p, sig2 = s2)


# Calculate the EWMA of Phase II observations.

n.phase2 <- 10000
set.seed(314)
X.phase2 <- matrix(rnorm(p * n.phase2), ncol = p, byrow = TRUE) %*% SO.rt

la <- 0.2
U.phase2 <- t(as.matrix(filter(la * X.phase2, 1 - la, method = "recursive")))

W.phase2 <- function(U, t, la, S.inv, S.irt, q) {
  stopifnot(la <= 1 && la > 0)
  stopifnot(is.integer(t))
  return((2 - la)/(la * (1 - (1 - la)^{2*t})) * W.phase1(X = U, S.inv = S.inv, S.irt = S.irt, q = q))
  #return(((2 - la)/la) * W.phase1(X = U, S.inv = S.inv, S.irt = S.irt, q = q))
}

W.phase2(U = U.phase2[, 10], t = as.integer(10), la = 0.2, S.inv = SO.inv, S.irt = SO.irt, q = 5)

ALO_phase2 <- function(t, U_prev, X_now, sig2, IC.bar, IC.sd, la, q) {
  stopifnot(length(U_prev) ==  length(X_now))
  stopifnot(length(IC.bar) == length(IC.sd))
  stopifnot((la > 0) && (la < 1))
  p <- length(U_prev)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  out <- .C("ALO_phase2", t = as.integer(t), p = p, U_prev = as.double(U_prev), X_now = as.double(X_now), 
            sig2 = as.double(sig2), IC.bar = as.double(IC.bar), IC.sd = as.double(IC.sd), la = as.double(la),
            U_now = double(p), W = double(p), q = q, Q = double(1))
  return(list(W = out$W, U = out$U_now, Q = out$Q))
}

ALO_phase2(t = 10, U_prev = U.phase2[, 9], X_now = X.phase2[10, ], sig2 = s2, IC.bar = ICstat[, "mean"],
           IC.sd = ICstat[, "std"], la = 0.2, q = 5)$W
ALO_phase2(t = 10, U_prev = U.phase2[, 9], X_now = X.phase2[10, ], sig2 = s2, IC.bar = ICstat[, "mean"],
           IC.sd = ICstat[, "std"], la = 0.2, q = 5)$U
U.phase2[, 10]


CC.W <- mcmapply(FUN = W.phase2, U = as.data.frame(U.phase2), t = seq_len(ncol(U.phase2)), 
                 MoreArgs = list(la = 0.2, S.inv = SO.inv, S.irt = SO.irt, q = 5), mc.cores = detectCores())
dim(CC.W)  

CC.Q <- function(W, IC.bar, IC.sd) {
  stopifnot(nrow(W) == length(IC.bar))
  stopifnot(length(IC.bar) == length(IC.sd))
  W.centered <- sweep(x = W, MARGIN = 1, STATS = IC.bar, FUN = "-")
  W.sdz <- sweep(x = W.centered, MARGIN = 1, STATS = IC.sd, FUN = "/")
  return(W.sdz)
}

W.sdz <- CC.Q(W = CC.W, IC.bar = ICstat[1:5, "mean"], IC.sd = ICstat[1:5, "std"])
QQ <- do.call(pmax.int, as.data.frame(t(W.sdz)))

QQ[3]
ALO_phase2(t = 3, U_prev = U.phase2[, 2], X_now = X.phase2[3, ], sig2 = s2, IC.bar = ICstat[, "mean"],
           IC.sd = ICstat[, "std"], la = 0.2, q = 5)$Q

RL <- function(Q, L) {
  return(min(c(which(Q > L)[1], length(Q)), na.rm = TRUE))
}

RL(Q = QQ, L = 5.181)

ALO_RL0 <- function(n.phase2, p, sig2, IC.bar, IC.sd, la, q, L){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("ALO_RL0", n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), IC.bar = as.double(IC.bar), 
            IC.sd = as.double(IC.sd), la = as.double(la), q = q, L = as.double(L), RL0 = as.integer(0))
  return(out$RL0)
}

set.seed(314)
ALO_RL0(n.phase2 = n.phase2, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], la = la, q = 3, L = 5.181)

ALO_Qs <- function(n.phase2, p, sig2, IC.bar, IC.sd, la, q){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("ALO_Qs", n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), IC.bar = as.double(IC.bar), 
            IC.sd = as.double(IC.sd), la = as.double(la), q = q, Qs = double(n.phase2))
  return(out$Qs)
}

set.seed(314)
alo.Qs <- ALO_Qs(n.phase2 = n.phase2, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], la = la, q = 3)
head(alo.Qs)

ALO_ARL0 <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q, L){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("ALO_ARL0", nrep = as.integer(nrep), n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), 
            IC.bar = as.double(IC.bar), IC.sd = as.double(IC.sd), la = as.double(la), q = q, L = as.double(L), 
            ARL0 = as.double(0.0))
  return(out$ARL0)
}

set.seed(314)
system.time(
  arl0 <- ALO_ARL0(nrep = 10000, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], 
                   la = la, q = p, L = 3.0)
  
)

ALO_Qmat <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("ALO_Qmat", nrep = as.integer(nrep), n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), 
            IC.bar = as.double(IC.bar), IC.sd = as.double(IC.sd), la = as.double(la), q = q, Qmat = double(nrep * n.phase2))
  return(out$Qmat)
}

set.seed(314)
system.time(
  alo.Qmat <- ALO_Qmat(nrep = 10000, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], 
                   la = la, q = p)
  
)

set.seed(22)
ALO_ARL0(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], la = la, q = p, L = 3.0)
set.seed(22)
alo.Qrep <- ALO_Qmat(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"], la = la, q = p)
head(alo.Qrep)
which.max(alo.Qrep > 3.0)

alo.Qmat <- matrix(alo.Qmat, nrow = 10000, byrow = TRUE)
alo.RLs <- apply(alo.Qmat > 3.0, 1, which.max)
mean(alo.RLs)
mm <- matrix(1:6, ncol = 3, byrow = TRUE)
apply(mm > 2.5, 1, which.max)


SCAD_ARL0 <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q, L, a = 3.7){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("SCAD_ARL0", nrep = as.integer(nrep), n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), a = as.double(a),
            IC.bar = as.double(IC.bar), IC.sd = as.double(IC.sd), la = as.double(la), q = q, L = as.double(L), ARL0 = as.double(0.0))
  return(out$ARL0)
}

set.seed(314)
system.time(
  scad_arl0 <- SCAD_ARL0(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = scadICstat[, "mean"], IC.sd = scadICstat[, "std"], 
                   la = la, q = p, L = 3.0)
  
)

SCAD_Qs <- function(n.phase2, p, sig2, IC.bar, IC.sd, la, q, a = 3.7){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("SCAD_Qs", n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), a = as.double(a), IC.bar = as.double(IC.bar), 
            IC.sd = as.double(IC.sd), la = as.double(la), q = q, Qs = double(n.phase2))
  return(out$Qs)
}

set.seed(314)
scad.Qs <- SCAD_Qs(n.phase2 = 2000, p = p, sig2 = s2, IC.bar = scadICstat[, "mean"], IC.sd = scadICstat[, "std"], la = la, q = p)
head(scad.Qs)
set.seed(314)
scad_arl0 <- SCAD_ARL0(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = scadICstat[, "mean"], IC.sd = scadICstat[, "std"], 
                       la = la, q = p, L = 3.0)

SCAD_Qmat <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q, a = 3.7){
  stopifnot(length(IC.bar) == length(IC.sd))
  p <- as.integer(p)
  stopifnot(p == length(IC.bar))
  q <- as.integer(q)
  stopifnot(q <= p)
  stopifnot(la > 0 & la <= 1)
  out <- .C("SCAD_Qmat", nrep = as.integer(nrep), n.phase2 = as.integer(n.phase2), p = p, sig2 = as.double(sig2), a = as.double(a),
            IC.bar = as.double(IC.bar), IC.sd = as.double(IC.sd), la = as.double(la), q = q, Qmat = double(nrep * n.phase2))
  return(out$Qmat)
}

set.seed(2200)
SCAD_ARL0(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = scadICstat[, "mean"], IC.sd = scadICstat[, "std"], la = la, q = p, L = 3.0)
set.seed(2200)
scad.Qrep <- SCAD_Qmat(nrep = 1, n.phase2 = 2000, p = p, sig2 = s2, IC.bar = scadICstat[, "mean"], IC.sd = scadICstat[, "std"], la = la, q = p)
head(scad.Qrep)


# search control limit (CL)

ALO_CL <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q, ARL0, lower, upper, ncpus = 4, precision = 0.2) {
  ncore <- min(ncpus, detectCores())
  stopifnot((nrep %% ncore) == 0)
  nrpercore <- as.integer(nrep / ncore)
  RLs <- mcmapply(ALO_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = upper),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  stopifnot(ARL > ARL0)
  RLs <- mcmapply(ALO_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = lower),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  stopifnot(ARL < ARL0)
  mid <- 0.5 * (lower + upper)
  RLs <- mcmapply(ALO_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = mid),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  while (abs(ARL - ARL0) >= precision) {
    if (ARL < ARL0) {
      lower <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(ALO_ARL0, nrep = rep(nrpercore, ncore), 
                      MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, 
                                      L = mid),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    } else {
      upper <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(ALO_ARL0, nrep = rep(nrpercore, ncore), 
                      MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, 
                                      L = mid),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    }
  }
  return(list(UCL = mid, ARL = ARL, sdARL = sd(RLs)/sqrt(ncore)))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  ALO_UCL.out <- ALO_CL(nrep = nrep, n.phase2 = 5000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"],
                        la = la, q = 3, ARL0 = 200, lower = 3, upper = 5, ncpus = 10, precision = 0.2)
)
invisible(rnorm(100)) # reset seed

SCAD_CL <- function(nrep, n.phase2, p, sig2, IC.bar, IC.sd, la, q, ARL0, lower, upper, a = 3.7, ncpus = 4, precision = 0.2) {
  ncore <- min(ncpus, detectCores())
  stopifnot((nrep %% ncore) == 0)
  nrpercore <- as.integer(nrep / ncore)
  RLs <- mcmapply(SCAD_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = upper,
                                  a = a),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  stopifnot(ARL > ARL0)
  RLs <- mcmapply(SCAD_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = lower, 
                                  a = a),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  stopifnot(ARL < ARL0)
  mid <- 0.5 * (lower + upper)
  RLs <- mcmapply(SCAD_ARL0, nrep = rep(nrpercore, ncore), 
                  MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, L = mid,
                                  a = a),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  while (abs(ARL - ARL0) >= precision) {
    if (ARL < ARL0) {
      lower <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(SCAD_ARL0, nrep = rep(nrpercore, ncore), 
                      MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, 
                                      L = mid, a = a),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    } else {
      upper <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(SCAD_ARL0, nrep = rep(nrpercore, ncore), 
                      MoreArgs = list(n.phase2 = n.phase2, p = p, sig2 = sig2, IC.bar = IC.bar, IC.sd = IC.sd, la = la, q = q, 
                                      L = mid, a = a),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    }
  }
  return(list(UCL = mid, ARL = ARL, sdARL = sd(RLs)/sqrt(ncore)))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  SCAD_UCL.out <- SCAD_CL(nrep = nrep, n.phase2 = 5000, p = p, sig2 = s2, IC.bar = ICstat[, "mean"], IC.sd = ICstat[, "std"],
                        la = la, q = 3, ARL0 = 200, lower = 3, upper = 5, a = 3.7, ncpus = 10, precision = 0.2)
)
invisible(rnorm(100)) # reset seed



# Wavelet representation

dyn.load("imspc.so")

wav2d_coef <- function(Z){
  n1 <- nrow(Z)
  n2 <- ncol(Z)
  J1 <- as.integer(log(n1, base = 2))
  J2 <- as.integer(log(n2, base = 2))
  p <- as.integer(2^{J1 + J2})
  out <- .C("wav2d_coef", n1 = n1, n2 = n2, Z = as.double(t(Z)), J1 = J1, J2 = J2, wcoef = double(p))
  return(list(J1 = out$J1, J2 = out$J2, wcoef = out$wcoef))
}

wav2d_img <- function(n1, n2, J1, J2, wcoef, thresh = 0) {
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  stopifnot(length(wcoef) == 2^{J1+J2})
  out <- .C("wav2d_img", n1 = as.integer(n1), n2 = as.integer(n2), J1 = as.integer(J1), J2 = as.integer(J2), 
            wcoef = as.double(wcoef), thresh = as.double(thresh), img = double(n1 * n2))
  x <- 0:(n1 - 1)/n1
  y <- 0:(n2 - 1)/n2
  img <- matrix(out$img, nrow = n1, byrow = TRUE)
  return(list(x = x, y = y, img = img))
}

J1 <- 7
J2 <- 7
wcoef <- 4^7:1
img0 <- wav2d_img(n1 = 2^J1, n2 = 2^J2, J1 = J1, J2 = J2, wcoef = wcoef)$img
wav_tf <- wav2d_coef(Z = img0)
wav_tf$J2
wcoef1 <- wav_tf$wcoef
img1 <- wav2d_img(n1 = 2^J1, n2 = 2^J2, J1 = J1, J2 = J2, wcoef = wcoef1)$img
max(abs(wcoef1 - wcoef))
max(abs(img0 - img1))

# IC images

n <- 128
f <- function(r) ifelse( r < 0.25, 1, 0)
x <- (0:(n - 1))/n
y <- x;
xy <- as.matrix(expand.grid(x, y))
rr <- sqrt((xy[, 1] - 0.5)^2 + (xy[, 2] - 0.5)^2)
ff <- matrix(f(rr), ncol = n)

image(ff, col = gray(c(0:255)/255))

beta_true <- wav2d_coef(Z = ff)$wcoef
p <- length(beta_true)

set.seed(314)
s2 <- 0.01
Z <- ff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n)
image(Z, col = gray(c(0:255)/255))

(SNR <- 10 * log(var(c(ff))/s2, base = 10))

wav_tf <- wav2d_coef(Z = Z)
Zhat <- wav2d_img(n1 = length(x), n2 = length(y), J1 = wav_tf$J1, J2 = wav_tf$J2, wcoef = wav_tf$wcoef, thresh = 0.002)$img
image(Zhat, col = gray(c(0:255)/255))

mean((Z - Zhat)^2)

# IC mean and std for image W's/V's: takes about 90 minutes to run

nrep <- 25000
ncores <- 10
nrpc <- as.integer(nrep / ncores)
RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  alo.mcout <- mcmapply(ALO_ICstat, n_IC = rep(nrpc, ncores), MoreArgs = list(p = length(beta_true), sig2 = s2/(n*n)), mc.cores = ncores)
)
invisible(rnorm(100)) # reset the random seed

alo.bar_mat <- alo.mcout[1:p, ]
alo.S2_mat <- alo.mcout[(p+1): (2*p), ]^2
alo.ic.bar <- rowMeans(alo.bar_mat)
alo.ic.std <- sqrt(rowMeans(alo.S2_mat) + rowMeans(alo.bar_mat^2) - alo.ic.bar^2)

head(cbind(alo.ic.bar, alo.ic.std))
max(alo.ic.std/sqrt(nrep))

RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  scad.mcout <- mcmapply(SCAD_ICstat, n_IC = rep(nrpc, ncores), MoreArgs = list(p = length(beta_true), sig2 = s2/(n*n), a = 3.7), 
                         mc.cores = ncores)
)
invisible(rnorm(100)) # reset the random seed

scad.bar_mat <- scad.mcout[1:p, ]
scad.S2_mat <- scad.mcout[(p+1): (2*p), ]^2
scad.ic.bar <- rowMeans(scad.bar_mat)
scad.ic.std <- sqrt(rowMeans(scad.S2_mat) + rowMeans(scad.bar_mat^2) - scad.ic.bar^2)

head(cbind(scad.ic.bar, scad.ic.std))
max(scad.ic.std/sqrt(nrep))

ICstats_df <- data.frame(alo.ic.bar = alo.ic.bar, alo.ic.std = alo.ic.std, scad.ic.bar = scad.ic.bar, scad.ic.std = scad.ic.std)
head(ICstats_df)

write.csv(x = ICstats_df, file = "ICstats_df.csv", row.names = FALSE)

## Read in the IC stats for ALASSO and SCAD charts.

dyn.load("imspc.so")

ICstats_df <- read.csv(file = "ICstats_df.csv")
head(ICstats_df)

alo.ic.bar <- ICstats_df[, "alo.ic.bar"]
alo.ic.std <- ICstats_df[, "alo.ic.std"]
scad.ic.bar <- ICstats_df[, "scad.ic.bar"]
scad.ic.std <- ICstats_df[, "scad.ic.std"]

# Phase II image monitoring

ALO_im <- function(t, U_prev, Z_now, wcoef_mean, sig2, IC.bar, IC.std, la, q) {
  p <- length(U_prev)
  stopifnot(p == length(wcoef_mean))
  stopifnot(p == length(IC.bar))
  stopifnot(p == length(IC.std))
  n1 <- nrow(Z_now)
  n2 <- ncol(Z_now)
  stopifnot(sig2 > 0)
  stopifnot(la > 0 && la <= 1)
  q <- as.integer(q)
  stopifnot(q <= p)
  out <- .C("ALO_im", t = as.integer(t), n1 = n1, n2 = n2, U_prev = as.double(U_prev), Z_now = as.double(t(Z_now)), 
            wcoef_mean = as.double(wcoef_mean), sig2 = as.double(sig2), IC.bar = as.double(IC.bar), 
            IC.std = as.double(IC.std), la = as.double(la), q = q, U_now = double(p), Q = double(1))
  return(list(U_now = out$U_now, Q = out$Q))
}

imm <- ALO_im(t = 1, U_prev = double(p), Z_now = Z, wcoef_mean = beta_true, sig2 = s2/(n*n), IC.bar = alo.ic.bar, IC.std = alo.ic.std,
              la = 0.2, q = p)
head(imm$U_now)
imm$Q

SCAD_im <- function(t, U_prev, Z_now, wcoef_mean, sig2, IC.bar, IC.std, la, q, a = 3.7) {
  p <- length(U_prev)
  stopifnot(p == length(wcoef_mean))
  stopifnot(p == length(IC.bar))
  stopifnot(p == length(IC.std))
  n1 <- nrow(Z_now)
  n2 <- ncol(Z_now)
  stopifnot(sig2 > 0)
  stopifnot(la > 0 && la <= 1)
  q <- as.integer(q)
  stopifnot(q <= p)
  out <- .C("SCAD_im", t = as.integer(t), n1 = n1, n2 = n2, U_prev = as.double(U_prev), Z_now = as.double(t(Z_now)), 
            wcoef_mean = as.double(wcoef_mean), sig2 = as.double(sig2), a = as.double(a), IC.bar = as.double(IC.bar), 
            IC.std = as.double(IC.std), la = as.double(la), q = q, U_now = double(p), Q = double(1))
  return(list(U_now = out$U_now, Q = out$Q))
}

SCAD_im(t = 1, U_prev = double(p), Z_now = Z, wcoef_mean = beta_true, sig2 = s2/(n*n), IC.bar = scad.ic.bar, IC.std = scad.ic.std,
       la = 0.2, q = p, a = 3.7)$Q

ALO_ims <- function(imseq, wcoef_mean, sig2, IC.bar, IC.std, la, q){
  p <- length(wcoef_mean)
  n_im <- dim(imseq)[3]
  Q <- double(n_im)
  U_prev <- double(p)
  for (t in 1:n_im) {
    mm <- ALO_im(t = t, U_prev = U_prev, Z_now = imseq[, , t], wcoef_mean = wcoef_mean, sig2 = sig2, IC.bar = IC.bar, 
                 IC.std = IC.std, la = la, q = q)
    Q[t] <- mm$Q
    U_prev <- mm$U_now
  }
  return(Q)
}

dyn.load("repsim.so")
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


ims <- array(NA, dim = c(n, n, 20))
fs <- 0.2
dff <- fault_square(n1 = n, n2 = n, size = fs)
#d.wcoef <- rep(0, length(beta_true))
#d.wcoef[c(2000, 10000)] <- 0.5
set.seed(314)
for (t in 1:20) {
  if (t < 15) {
    ims[, , t] <- ff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n, byrow = TRUE)  
  } else {
    ims[, , t] <- ff + dff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n, byrow = TRUE)  
    # ims[, , t] <- ff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n)
    # ims[50:55, 70:75 , t] <- ims[50:55, 70:75, t] + 0.5
    # ims[20:25, 16:20 , t] <- ims[20:25, 16:20, t] - 0.4
    # ims[, , t] <- wav2d_img(n1 = as.integer(n), n2 = as.integer(n), J1 = as.integer(log(n, base = 2)), 
    #                         J2 = as.integer(log(n, base = 2)), wcoef = as.double(beta_true + d.wcoef), 
    #                         thresh = as.double(0))$img + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n)
  }
}

layout(matrix(c(1:6, rep(7, 6)), 2, 6, byrow = TRUE))
for (t in 12: 17) {
  image(ims[, , t], col = gray(c(0:255)/255), main = paste0("Image #", t), zlim = range(ims))
}

## make figure 1 for publication

par(mfrow = c(2, 3), mar = c(2, 3, 1, 3))
image(ff, col = gray.colors(256), xaxt = "n", yaxt = "n", zlim = range(Z))
title(sub = "(a)", line = 0.2, cex = 2)
image(Z, col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(b)", line = 0.2, cex = 2)
image(Z + fault_cross(n1 = n, n2 = n, size = 0.6), col = gray.colors(256), xaxt = "n", yaxt = "n") # alter the size to make the fault visible.
title(sub = "(c)", line = 0.2, cex = 2)
image(Z + fault_cross(n1 = n, n2 = n, size = 1.0), col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(d)", line = 0.2, cex = 2)
image(Z + fault_square(n1 = n, n2 = n, size = 0.2), col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(e)", line = 0.2, cex = 2)
image(Z + fault_square(n1 = n, n2 = n, size = 0.4), col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(f)", line = 0.2, cex = 2)

## end of making figure 1


alo.Qs <- ALO_ims(imseq = ims, wcoef_mean = beta_true, sig2 = s2/(n*n), IC.bar = alo.ic.bar, IC.std = alo.ic.std, la = 0.2, q = p)
par(mfrow = c(1, 1))
plot(alo.Qs, type = "b")
alo.Qs[1:20]

SCAD_ims <- function(imseq, wcoef_mean, sig2, IC.bar, IC.std, la, q, a = 3.7){
  p <- length(wcoef_mean)
  n_im <- dim(imseq)[3]
  Q <- double(n_im)
  U_prev <- double(p)
  for (t in 1:n_im) {
    mm <- SCAD_im(t = t, U_prev = U_prev, Z_now = imseq[, , t], wcoef_mean = wcoef_mean, sig2 = sig2, IC.bar = IC.bar, 
                 IC.std = IC.std, la = la, q = q, a = a)
    Q[t] <- mm$Q
    U_prev <- mm$U_now
  }
  return(Q)
}

layout(matrix(c(1:6, rep(7, 6)), 2, 6, byrow = TRUE))
for (t in 12: 17) {
  image(ims[, , t], col = gray(c(0:255)/255), main = paste0("Image #", t), zlim = range(ims))
}

scad.Qs <- SCAD_ims(imseq = ims, wcoef_mean = beta_true, sig2 = s2/(n*n), IC.bar = scad.ic.bar, IC.std = scad.ic.std, la = 0.2, q = p, a = 3.7)
plot(scad.Qs, type = "b")
scad.Qs[1:20]

## Post signal diagnostics

change_loc <- function(ims, wcoef_mean) {
  n_im <- dim(ims)[3]
  n1 <- dim(ims)[1]
  n2 <- dim(ims)[2]
  p <- length(wcoef_mean)
  stopifnot(p == as.integer(2^{as.integer(log(n1 * n2, base = 2))}))
  Z_seq <- as.double(c(aperm(ims, perm = c(2, 1, 3))))
  out <- .C("change_loc", n_im = n_im, n1 = n1, n2 = n2, Z_seq = Z_seq, wcoef_mean = as.double(wcoef_mean), 
            like = double(n_im), loc_hat = as.integer(0))
  return(list(like = out$like, loc = out$loc_hat))
}

system.time(print(change_loc(ims = ims[, , 1:17], wcoef_mean = beta_true)))

ALO_change_wcoef <- function(ims.OC, sig2, wcoef_mean){
  stopifnot(sig2 > 0)
  n_im <- dim(ims.OC)[3]
  n1 <- dim(ims.OC)[1]
  n2 <- dim(ims.OC)[2]
  p <- length(wcoef_mean)
  stopifnot(p == as.integer(2^{as.integer(log(n1 * n2, base = 2))}))
  Z_seq <- as.double(c(aperm(ims.OC, perm = c(2, 1, 3))))
  out <- .C("ALO_change_wcoef", n_im = n_im, n1 = n1, n2 = n2, sig2 = as.double(sig2), Z_seq = Z_seq, 
            wcoef_mean = as.double(wcoef_mean), EBIC = double(p), OC_wcoef = double(p))
  return(list(EBIC = out$EBIC, OC_wcoef = out$OC_wcoef))
}

SCAD_change_wcoef <- function(ims.OC, sig2, wcoef_mean, a = 3.7){
  stopifnot(sig2 > 0)
  n_im <- dim(ims.OC)[3]
  n1 <- dim(ims.OC)[1]
  n2 <- dim(ims.OC)[2]
  p <- length(wcoef_mean)
  stopifnot(p == as.integer(2^{as.integer(log(n1 * n2, base = 2))}))
  Z_seq <- as.double(c(aperm(ims.OC, perm = c(2, 1, 3))))
  out <- .C("SCAD_change_wcoef", n_im = n_im, n1 = n1, n2 = n2, sig2 = as.double(sig2), a = as.double(a), Z_seq = Z_seq, 
            wcoef_mean = as.double(wcoef_mean), EBIC = double(p), OC_wcoef = double(p))
  return(list(EBIC = out$EBIC, OC_wcoef = out$OC_wcoef))
}

alo.OC.out <- ALO_change_wcoef(ims.OC = ims[, , 15:17], sig2 = s2/(n*n), wcoef_mean = beta_true)
which.min(alo.OC.out$EBIC)
which(abs(alo.OC.out$OC_wcoef) > 0)
scad.OC.out <- SCAD_change_wcoef(ims.OC = ims[, , 15:17], sig2 = s2/(n*n), wcoef_mean = beta_true, a = 3.7)
which.min(scad.OC.out$EBIC)
which(abs(scad.OC.out$OC_wcoef) > 0)

print(alo.OC.out$EBIC[1:6], digits = 10)
print(scad.OC.out$EBIC[1:6], digits = 10)

J1 <- as.integer(log(n, base = 2))
J2 <- as.integer(log(n, base = 2))
img.diag <- wav2d_img(n1 = n, n2 = n, J1 = J1, J2 = J2, wcoef = alo.OC.out$OC_wcoef)$img
img.diag <- wav2d_img(n1 = n, n2 = n, J1 = J1, J2 = J2, wcoef = scad.OC.out$OC_wcoef)$img
img.trueOC <- wav2d_img(n1 = n, n2 = n, J1 = J1, J2 = J2, wcoef = beta_true + d.wcoef)$img

image(img.diag, col = gray(c(0:255)/255))

par(mfrow = c(2, 2))
image(ims[, , 14], col = gray(c(0:255)/255), zlim = range(c(c(ims), c(img.diag))))
image(ims[, , 15], col = gray(c(0:255)/255), zlim = range(c(c(ims), c(img.diag))))
image(abs(img.diag) > 0, col = gray(c(0:255)/255))
image(dff, col = gray(c(0:255)/255))

OC.out$OC_wcoef[which(abs(OC.out$OC_wcoef) > 0)]
d.wcoef[which(abs(d.wcoef) > 0)]
sum(abs(OC.out$OC_wcoef - d.wcoef))
sum(abs(img.diag - img.trueOC))
range(img.diag)
range(ims)
range(img.trueOC)

## From 2-d wavelet coefficients to pixel locations

dyn.load("imspc.so")
n <- 128
d.wcoef <- rep(0, n^2)
d.wcoef[c(2000, 10000)] <- 0.5
J1 <- 7
J2 <- 7
img0 <- wav2d_img(n1 = 2^J1, n2 = 2^J2, J1 = J1, J2 = J2, wcoef = d.wcoef)$img
img0 <- ifelse(abs(img0) > 0, 1, 0)
image(img0, col = gray(c(0:255)/255))

which(img0 > 0, arr.ind = TRUE)

w2d_pxl <- function(n1, n2, J1, J2, w2d_indx) {
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  J1 <- as.integer(J1)
  J2 <- as.integer(J2)
  stopifnot(2^J1 <= n1)
  stopifnot(2^J2 <= n2)
  p <- 2^{J1 + J2}
  stopifnot(max(w2d_indx) <= p)
  pxl <- matrix(.C("w2d_pxl", n1 = n1, n2 = n2, J1 = J1, J2 = J2, nb = length(w2d_indx), bseq = as.integer(w2d_indx - 1), 
                   pxl = integer(n1 * n2))$pxl, nrow = n1, byrow = TRUE)
  return(pxl)
}

pxl <- w2d_pxl(n1 = n, n2 = n, J1 = J1, J2 = J2, w2d_indx = which(abs(d.wcoef) > 0))
sum(abs(pxl - img0))

ic.img <- matrix(0, nrow = n, ncol = n)
oc.img <- ic.img
oc.img[10:11, 17:32] <- 1
diff.img <- oc.img - ic.img
image(ic.img, col = gray(c(0:255)/255))
image(oc.img, col = gray(c(0:255)/255))
diff.wcoef <- wav2d_coef(Z = diff.img)$wcoef
oc.pxl <- wav2d_img(n1 = n, n2 = n, J1 = J1, J2 = J2, wcoef = diff.wcoef)$img
oc.pxl[abs(oc.pxl) > 1/10^{15}] <- 1
oc.pxl[abs(oc.pxl) < 1] <- 0
image(oc.pxl, col = gray(c(0:255)/255))
sum(abs(diff.img - oc.pxl))

## dQ measure

nn <- 128
Omega <- matrix(0, nrow = nn, ncol = nn)
D <- Omega
Dhat <- Omega
D[35:37, 1:10] <- 1
Dhat[35:37, 6:15] <- 1
w <- 0.5
(w * sum((Dhat - D) > 0)/(nn^2 - sum(D)) + (1 - w) * sum((D - Dhat) > 0) / sum(D))
.C("dQ", n1 = as.integer(nn), n2 = as.integer(nn), D = as.integer(t(D)), Dhat = as.integer(t(Dhat)), w = as.double(w), dq = double(1))$dq


## stacked 1-d wavelet representations of images

dyn.load("imspc.so")
n <- 128
zz <- matrix(0, nrow = n, ncol = n)
zz[10:11, 17:32] <- 1
image(zz, col = gray(c(0:255)/255))

wav1ds <- function(Z, decom.level = 4){
  n1 <- nrow(Z)
  n2 <- ncol(Z)
  J <- as.integer(decom.level)
  stopifnot(2^J <= n2)
  Phi <- t(matrix(.C("wav1d_dm", n = n2, J = J, dm = double(2^J * n2))$dm, ncol = n2, byrow = TRUE))
  BETA <- Z %*% Phi / n2
  Zhat <- BETA %*% t(Phi)
  return(list(BETA = BETA, Zhat = Zhat))
}

zz.hat <- wav1ds(Z = zz, decom.level = 7)$Zhat
image(zz.hat, col = gray(c(0:255)/255))
sum(abs(zz - zz.hat))

zz.BETA <- wav1ds(Z = zz, decom.level = 7)$BETA

wav1ds_img <- function(BETA, n2) {
  n1 <- nrow(BETA)
  n2 <- as.integer(n2)
  stopifnot(ncol(BETA) <= n2)
  J <- as.integer(log(ncol(BETA), base = 2))
  Phi <- t(matrix(.C("wav1d_dm", n = n2, J = J, dm = double(2^J * n2))$dm, ncol = n2, byrow = TRUE))
  img <- BETA %*% t(Phi)
  return(img)
}

zz.img <- wav1ds_img(BETA = zz.BETA, n2 = ncol(zz))
image(zz.img, col = gray(c(0:255)/255))
sum(abs(zz - zz.img))

dBETA <- matrix(0, nrow = nrow(zz.BETA), ncol = ncol(zz.BETA))
dBETA[100:101, 33:48] <- 1
image(dBETA, col = gray(c(0:255)/255))
oc.img <- wav1ds_img(BETA = zz.BETA + dBETA, n2 = ncol(zz))
image(oc.img, col = gray(c(0:255)/255))
image(abs(oc.img) > 0, col = gray.colors(256))

## mapping between stacked 1-d wavelet coefficients and pixels.

wav1ds_pxl <- function(n1, n2, ocBETA, eps = 1/10^{15}) {
  stopifnot(n1 == nrow(ocBETA))
  n1 <- nrow(ocBETA)
  n2 <- as.integer(n2)
  stopifnot(ncol(ocBETA) <= n2)
  J <- as.integer(round(log(ncol(ocBETA), base = 2)))
  ocBETA01 <- matrix(0, nrow = n1, ncol = ncol(ocBETA))
  ocBETA01[abs(ocBETA) > eps] <- 1
  pxl <- .C("w1ds_pxl", n1 = n1, n2 = n2, J2 = J, BETA01 = as.integer(t(ocBETA01)), pxl = integer(n1 * n2))$pxl
  pxl <- matrix(pxl, nrow = n1, byrow = TRUE)
  return(pxl)
}

oc.pxl <- wav1ds_pxl(n1 = nrow(zz.img), n2 = ncol(zz.img), ocBETA = dBETA)
image(oc.pxl, col = gray(c(0:255)/255))
diff.img <- oc.img - zz.img
diff.img[abs(diff.img) > 0] <- 1
diff.img[abs(diff.img) < 1] <- 0
image(diff.img, col = gray(c(0:255)/255))
sum(abs(diff.img - oc.pxl))

## check: detecting OC pixels

ic.img <- matrix(0, nrow = n, ncol = n)
oc.img <- ic.img
oc.img[10:11, 17:32] <- 1
diff.img <- oc.img - ic.img
image(ic.img, col = gray(c(0:255)/255))
image(oc.img, col = gray(c(0:255)/255))
ic.BETA <- wav1ds(Z = ic.img, decom.level = 7)$BETA
oc.BETA <- wav1ds(Z = oc.img, decom.level = 7)$BETA
dBETA <- oc.BETA - ic.BETA
d.img <- wav1ds_img(BETA = dBETA, n2 = ncol(diff.img))
oc.pxl <- d.img
oc.pxl[abs(oc.pxl) > 1/10^{15}] <- 1
oc.pxl[abs(oc.pxl) < 1] <- 0
image(oc.pxl, col = gray(c(0:255)/255))
sum(abs(oc.pxl - diff.img))

## GLR IC control limit

dyn.load("imspc.so")

GLR1d_RLs0 <- function(nRuns, upbnd, wsize, UCL) {
  nRuns <- as.integer(nRuns)
  upbnd <- as.integer(upbnd)
  m <- as.integer(wsize)
  UCL <- as.double(UCL)
  out <- .C("GLR1d_RLs0", nRuns = nRuns, upbnd = upbnd, m = m, UCL = UCL, RLs0 = double(nRuns))
  return(out$RLs0)
}

set.seed(314)
system.time(RL <- GLR1d_RLs0(nRuns = 10^4, upbnd = 2 * 10^4, wsize = 10, UCL = 5.0669)) 
mean(RL)
sd(RL)/sqrt(10^4)

IdpGLR_RLs0 <- function(nGLRs, nRuns, upbnd, wsize, UCL) {
  nGLRs <- as.integer(nGLRs)
  nRuns <- as.integer(nRuns)
  upbnd <- as.integer(upbnd)
  m <- as.integer(wsize)
  UCL <- as.double(UCL)
  out <- .C("IdpGLR_RLs0", nGLRs = nGLRs, nRuns = nRuns, upbnd = upbnd, m = m, UCL = UCL, RLs0 = double(nRuns))
  return(out$RLs0)
}

set.seed(314)
IdpRL <- IdpGLR_RLs0(nGLRs = 1, nRuns = 10^4, upbnd = 2 * 10^4, wsize = 10, UCL = 5.0669) 
mean(IdpRL)

set.seed(314)
system.time(IdpRL <- IdpGLR_RLs0(nGLRs = 128 * 2^4, nRuns = 10, upbnd = 200 * 10, wsize = 10, UCL = 12.75))
mean(IdpRL)
sd(IdpRL)/sqrt(10)

m <- 10
nrpercore <- 1000
ncore <- 12
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(314) # takes 5 minutes to run
system.time(
  IdpRLs <- mcmapply(IdpGLR_RLs0, nRuns = rep(nrpercore, ncore), MoreArgs = list(nGLRs = 128 * 2^4, upbnd = 200 * 10, wsize = m,
                                                                                 UCL = 12.75),
                     mc.cores = detectCores())
)
invisible(rnorm(100)) # reset seed
mean(IdpRLs)
sd(IdpRLs)/sqrt(nrpercore * ncore)

GLR_CL <- function(nGLRs, nRuns, upbnd, wsize, ARL0, lower, upper, ncpus = 4, precision = 0.2) {
  ncore <- min(ncpus, detectCores())
  nrpercore <- as.integer(nRuns/ncore)
  m <- as.integer(wsize)
  mid <- 0.5 * (lower + upper)
  RLs <- mcmapply(IdpGLR_RLs0, nRuns = rep(nrpercore, ncore), MoreArgs = list(nGLRs = nGLRs, upbnd = upbnd, wsize = m, UCL = mid),
                  mc.cores = ncore)
  ARL <- mean(RLs)
  while (abs(ARL - ARL0) >= precision) {
    if (ARL < ARL0) {
      lower <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(IdpGLR_RLs0, nRuns = rep(nrpercore, ncore), MoreArgs = list(nGLRs = nGLRs, upbnd = upbnd, wsize = m, UCL = mid),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    } else {
      upper <- mid
      mid <- 0.5 * (lower + upper)
      RLs <- mcmapply(IdpGLR_RLs0, nRuns = rep(nrpercore, ncore), MoreArgs = list(nGLRs = nGLRs, upbnd = upbnd, wsize = m, UCL = mid),
                      mc.cores = ncore)
      ARL <- mean(RLs)
    }
  }
  return(list(UCL = mid, ARL = ARL, sdARL = sd(RLs)/sqrt(nRuns)))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(42)
system.time(
  UCL.out <- GLR_CL(nGLRs = 1, nRuns = 10^5, upbnd = 2 * 10^4, wsize = 10, ARL0 = 200, lower = 4, upper = 6, ncpus = 10)
)
invisible(rnorm(100)) # reset seed

mean(IdpGLR_RLs0(nGLRs = 1, nRuns = 10^5, upbnd = 2 * 10^4, wsize = 10, UCL = UCL.out$UCL))

## GLR monitoring

dyn.load("imspc.so")

nGLRs <- 5
m <- 10
sigma <- 0.1
set.seed(42)
Xmat <- matrix(rnorm(nGLRs * m, sd = sigma), ncol = m)
#Xmat[, 5:m] <- 0
BS <- t(apply(Xmat[, m:1], 1, FUN = cumsum))[, m:1]
set.seed(737)
X <- rnorm(nGLRs, sd = sigma)
Xnew <- as.matrix(cbind(Xmat, X))
Xnew <- t(apply(Xnew[, (m + 1) : 1], 1, FUN = cumsum))[, (m + 1): 1]
#Xnew <- Xmat
#Xnew[, 5] <- X
#BSnew <- t(apply(Xnew[, (m) : 1], 1, FUN = cumsum))[, (m): 1]
BSnew <- Xnew[, -1]
GLRmat <- BSnew^2 %*% diag(1/c(m:1)) / (sigma^2) * 0.5
#GLRmat <- BSnew^2 %*% diag(1/c(5:1, rep(1, 5))) / (sigma^2) * 0.5
GLRstats <- apply(GLRmat, 1, FUN = max)
GLRtaus <- apply(GLRmat, 1, FUN = which.max)
s <- 11
#s <- 5
out <- .C("GLRs_phase2", nGLRs = as.integer(nGLRs), m = as.integer(m), s = as.integer(s), sig2 = as.double(sigma^2), 
          BSprev = as.double(t(BS)), Xnow = as.double(X), Xmean = double(nGLRs), BSnow = double(nGLRs * m), GLRstats = double(nGLRs), 
          GLRtaus = integer(nGLRs))
out$GLRstats
s - (m - (out$GLRtaus + 1))

GLR_im <- function(t, BSprev, Xnow, Xmean, betasig2) {
  stopifnot(t >= 0)
  stopifnot(betasig2 > 0)
  nGLRs <- nrow(BSprev)
  m <- ncol(BSprev)
  stopifnot(nGLRs == length(Xnow))
  stopifnot(nGLRs == length(Xmean))
  out <- .C("GLRs_phase2", nGLRs = nGLRs, m = m, s = as.integer(t), sig2 = as.double(betasig2), BSprev = as.double(t(BSprev)), Xnow = as.double(Xnow), 
            Xmean = as.double(Xmean), BSnow = double(nGLRs * m), GLRstats = double(nGLRs), GLRtaus = integer(nGLRs))
  if (t < m) {
    GLRtaus <- out$GLRtaus
  } else {
    GLRtaus = t - m + out$GLRtaus
  }
  BSnow <- matrix(out$BSnow, ncol = m, byrow = TRUE)
  return(list(BSnow = BSnow, GLRstats = out$GLRstats, GLRtaus = GLRtaus))
}

GLRim.out <- GLR_im(t = s, BSprev = BS, Xnow = X, Xmean = rep(0, length(X)), betasig2 = sigma^2)
GLRim.out$GLRstats
GLRim.out$GLRtaus

GLR_ims <- function(imseq, decom.level, betasig2, m, wcoef_mean) {
  nr <- dim(imseq)[1]
  nc <- dim(imseq)[2]
  n_im <- dim(imseq)[3]
  J <- as.integer(decom.level)
  stopifnot(2^J <= nc)
  nGLRs <- as.integer(nr * 2^J)
  m <- as.integer(m)
  stopifnot(nGLRs == length(wcoef_mean))
  GLRchart <- rep(NA, n_im)
  GLRstatsMat <- matrix(0, nrow = nGLRs, ncol = n_im)
  GLRtausMat <- matrix(0, nrow = nGLRs, ncol = n_im)
  BSprev <- matrix(0, nrow = nGLRs, ncol = m)
  for (t in 1:n_im) {
    betastack <- c(t(wav1ds(Z = imseq[, , t], decom.level = J)$BETA))
    GLRupdate <- GLR_im(t = t, BSprev = BSprev, Xnow = betastack, Xmean = wcoef_mean, betasig2 = betasig2)
    BSprev <- GLRupdate$BSnow
    GLRchart[t] <- max(GLRupdate$GLRstats)
    GLRstatsMat[, t] <- GLRupdate$GLRstats
    GLRtausMat[, t] <- GLRupdate$GLRtaus
  }
  return(list(GLRchart = GLRchart, GLRstatsMat = GLRstatsMat, GLRtausMat = GLRtausMat))
}

dyn.load("repsim.so")
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

n_ims <- 20
ims <- array(NA, dim = c(n, n, n_ims))
fs <- 0.4
dff <- fault_cross(n1 = n, n2 = n, size = fs)
dff <- fault_square(n1 = n, n2 = n, size = fs)
#d.wcoef <- rep(0, length(beta_true))
#d.wcoef[c(2000, 10000)] <- 0.5
set.seed(314)
for (t in seq_len(n_ims)) {
  if (t < 15) {
    ims[, , t] <- ff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n, byrow = TRUE)  
  } else {
    ims[, , t] <- ff + dff + matrix(rnorm(length(x) * length(y), sd = sqrt(s2)), ncol = n, byrow = TRUE)  
  }
}

layout(matrix(c(1:2, rep(3, 2)), 2, 2, byrow = TRUE))
for (t in 14: 15) {
  image(ims[, , t], col = gray(c(0:255)/255), main = paste0("Image #", t), zlim = range(ims))
}

J <- 4
betastack_mean <- c(t(wav1ds(Z = ff, decom.level = J)$BETA))
GLRs.out <- GLR_ims(imseq = ims, decom.level = J, betasig2 = s2/n, m = 10, wcoef_mean = betastack_mean)

plot(GLRs.out$GLRchart, type = "b")

all(GLRs.out$GLRtausMat[, 1] == 0)

which(GLRs.out$GLRchart >= 12.75)
GLRs.out$GLRtausMat[which.max(GLRs.out$GLRstatsMat[, 15]), 15]

Rs <- GLRs.out$GLRstatsMat[, 15]
Zoc1st <- ims[, , 15]
WToc1st <- wav1ds(Z = Zoc1st, decom.level = 4)
Boc1st <- c(t(WToc1st$BETA))
Boc1st[Rs < 12.75] <- 0
Boc1st <- matrix(Boc1st, nrow = nrow(Zoc1st), byrow = TRUE)
PXLoc1st <- wav1ds_img(BETA = Boc1st, n2 = ncol(Zoc1st))
image(PXLoc1st, col = gray(c(0:255)/255))

## Repeated Simulations

dyn.load("repsim.so")

imgen <- function(n1, n2) {
  out <- matrix(.C("imgen", n1 = as.integer(n1), n2 = as.integer(n2), ff = double(n1 * n2))$ff, ncol = n2, byrow = TRUE)
  return(out)
}

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

## set parameters

n <- 128
s2 <- 0.01
fs <- 0.4
tau <- 14
arl0 <- 200
rlbnd <- 10 * 200
NRep <- 10^4
la <- 0.2
p <- 2^{as.integer(log(n, base = 2)) * 2}
q <- p

img0 <- imgen(n1 = n, n2 = n)
image(img0, col = gray(c(0:255)/255))
cross1 <- fault_cross(n1 = n, n2 = n, size = fs)
square1 <- fault_square(n1 = n, n2 = n, size = fs)
set.seed(314)
noise <- matrix(rnorm(n * n, sd = sqrt(s2)), ncol = n)
par(mfrow = c(2, 2))
image(img0 + cross1, col = gray(c(0:255)/255))
image(img0 + square1, col = gray(c(0:255)/255))
image(img0 + cross1 + noise, col = gray(c(0:255)/255))
image(img0 + square1 + noise, col = gray(c(0:255)/255))

## Repeat ALASSO

### ALASSO control limit from the big run: 4.028931

set.seed(314)
system.time(
  print(ALO_ARL0(nrep = 1, n.phase2 = rlbnd, p = p, sig2 = s2/(n^2), IC.bar = alo.ic.bar, IC.sd = alo.ic.std, 
                   la = la, q = p, L = 4.028931))
  
)

library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  RLs <- mcmapply(ALO_ARL0, nrep = rep(1, 10), 
                  MoreArgs = list(n.phase2 = rlbnd, p = p, sig2 = s2/(n^2), IC.bar = alo.ic.bar, IC.sd = alo.ic.std, la = la, q = p, L = 4.028931),
                  mc.cores = 10)
  
)
invisible(rnorm(100)) # reset seed
mean(RLs)

Qmat.alo <- read.csv(file = "QMalo.csv", header = TRUE)
Qmat.alo <- as.matrix(Qmat.alo)
RLs.alo <- apply(Qmat.alo > 4.025, 2, FUN = which.max)
mean(RLs.alo)
sd(RLs.alo)/sqrt(ncol(Qmat.alo))

QsRLs <- function(QM, L) {
  QM01 <- (as.matrix(QM) > L)
  n1 <- nrow(QM)
  n2 <- ncol(QM)
  QM01[n1, ] <- rep(TRUE, n2)
  RLs <- apply(QM01, 2, FUN = which.max)
  return(RLs)
}
mean(QsRLs(QM = Qmat.alo, L = 4.025))

### ALASSO phase II simulations

ALO_rep <- function(nrep, n1, n2, faultshape, faultsize, noiselevel, tau, RLbnd, UCL, IC.bar, IC.std, la, q) {
  J1 <- as.integer(log(n1, base = 2))
  J2 <- as.integer(log(n2, base = 2))
  p <- 2^{J1 + J2}
  stopifnot(q <= p)
  stopifnot(p == length(IC.bar))
  stopifnot(p == length(IC.std))
  stopifnot(noiselevel > 0)
  stopifnot(la > 0 && la <= 1)
  if (faultshape == "cross") {
    shape <- as.integer(0)
  } else {
    if (faultshape == "square") {
      shape <- as.integer(1)
    } else {
      stop("fault shape has be either 'cross' or 'shape'. ")
    }
  }
  out <- .C("ALO_rep", nrep = as.integer(nrep), n1 = as.integer(n1), n2 = as.integer(n2), shape = shape, faultsize = as.double(faultsize), 
            noiselevel = as.double(noiselevel), tau = as.integer(tau), RLbnd = as.integer(RLbnd), UCL = as.double(UCL), 
            IC.bar = as.double(IC.bar), IC.std = as.double(IC.std), la = as.double(la), q = as.integer(q), RL1 = integer(nrep), 
            tauhat = integer(nrep), ocPXL = double(n1 * n2), dQs = double(nrep))
  ocPXL <- matrix(out$ocPXL, ncol = n2, byrow = TRUE)
  return(list(RL1 = out$RL1, tauhat = out$tauhat, ocPXL = ocPXL, dQs = out$dQs))
}

set.seed(314)
system.time(
  ALOrep.out <- ALO_rep(nrep = 1, n1 = n, n2 = n, faultshape = "cross", faultsize = fs+0*sqrt(s2), noiselevel = sqrt(s2), tau = tau, RLbnd = rlbnd, 
                      UCL = 4.1, IC.bar = alo.ic.bar, IC.std = alo.ic.std, la = 0.2, q = n^2)
)
mean(ALOrep.out$RL1)
mean(ALOrep.out$tauhat)
image(ALOrep.out$ocPXL, col = gray(c(0:255)/255))
Dhat <- matrix(0, nrow = n, ncol = n)
Dhat[which(abs(ALOrep.out$ocPXL) > 10^{-15})] <- 1
image(Dhat, col = gray(c(0:255)/255))
D <- matrix(0, nrow = n, ncol = n)
D[which(abs(cross1) > 10^{-15})] <- 1
image(D, col = gray(c(0:255)/255))
(0.5 * sum((Dhat - D) > 0)/(n^2 - sum(D)) + (1 - 0.5) * sum((D - Dhat) > 0) / sum(D))
mean(ALOrep.out$dQs)
max(abs(img.diag - ALOrep.out$ocPXL))

## Repeat SCAD

### SCAD control limit

set.seed(314)
system.time(
  print(SCAD_ARL0(nrep = 10, n.phase2 = rlbnd, p = p, sig2 = s2/(n^2), IC.bar = scad.ic.bar, IC.sd = scad.ic.std, 
                 la = la, q = p, L = 4.0))
  
)

Qmat.scad <- read.csv(file = "QMscad.csv", header = TRUE)
Qmat.scad <- as.matrix(Qmat.scad)
RLs.scad <- apply(Qmat.scad > 4.0, 2, FUN = which.max)
mean(RLs.scad)
sd(RLs.scad)/sqrt(ncol(Qmat.scad))

RLs.scad <- QsRLs(QM = Qmat.scad, L = 4.0)
c(mean(RLs.scad), sd(RLs.scad)/sqrt(ncol(Qmat.scad)))

### SCAD phase II simulations

SCAD_rep <- function(nrep, n1, n2, faultshape, faultsize, noiselevel, tau, RLbnd, UCL, IC.bar, IC.std, la, q, a = 3.7) {
  J1 <- as.integer(log(n1, base = 2))
  J2 <- as.integer(log(n2, base = 2))
  p <- 2^{J1 + J2}
  stopifnot(q <= p)
  stopifnot(p == length(IC.bar))
  stopifnot(p == length(IC.std))
  stopifnot(noiselevel > 0)
  stopifnot(la > 0 && la <= 1)
  if (faultshape == "cross") {
    shape <- as.integer(0)
  } else {
    if (faultshape == "square") {
      shape <- as.integer(1)
    } else {
      stop("fault shape has be either 'cross' or 'shape'. ")
    }
  }
  out <- .C("SCAD_rep", nrep = as.integer(nrep), n1 = as.integer(n1), n2 = as.integer(n2), shape = shape, faultsize = as.double(faultsize), 
            noiselevel = as.double(noiselevel), tau = as.integer(tau), RLbnd = as.integer(RLbnd), UCL = as.double(UCL), a = as.double(a), 
            IC.bar = as.double(IC.bar), IC.std = as.double(IC.std), la = as.double(la), q = as.integer(q), RL1 = integer(nrep), 
            tauhat = integer(nrep), ocPXL = double(n1 * n2), dQs = double(nrep))
  ocPXL <- matrix(out$ocPXL, ncol = n2, byrow = TRUE)
  return(list(RL1 = out$RL1, tauhat = out$tauhat, ocPXL = ocPXL, dQs = out$dQs))
}

set.seed(314)
system.time(
  SCADrep.out <- SCAD_rep(nrep = 1, n1 = n, n2 = n, faultshape = "cross", faultsize = fs+0*sqrt(s2), noiselevel = sqrt(s2), tau = tau, RLbnd = rlbnd, 
                        UCL = 4.0, IC.bar = scad.ic.bar, IC.std = scad.ic.std, la = 0.2, q = n^2)
)

mean(SCADrep.out$RL1)
mean(SCADrep.out$tauhat)
max(abs(img.diag - SCADrep.out$ocPXL))

## Repeat GLR

### GLR control limit

set.seed(314)
system.time(IdpRL <- IdpGLR_RLs0(nGLRs = 128 * 2^4, nRuns = 1000, upbnd = 200 * 10, wsize = 10, UCL = 12.75))
mean(IdpRL)
sd(IdpRL)/sqrt(1000)

### GLR phase II simulations
dyn.load("imspc.so")
dyn.load("repsim.so")

GLR_rep <- function(nrep, n1, n2, faultshape, faultsize, noiselevel, tau, RLbnd, UCL, J2, m = 10) {
  nrep <- as.integer(nrep)
  m <- as.integer(m)
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)
  J2 <- as.integer(J2)
  p1d <- as.integer(2^{J2})
  stopifnot(p1d <= n2)
  stopifnot(noiselevel > 0)
  if (faultshape == "cross") {
    shape <- as.integer(0)
  } else {
    if (faultshape == "square") {
      shape <- as.integer(1)
    } else {
      stop("fault shape has be either 'cross' or 'shape'. ")
    }
  }
  out <- .C("GLR_rep", nrep = nrep, n1 = n1, n2 = n2, shape = shape, faultsize = as.double(faultsize), 
            noiselevel = as.double(noiselevel), tau = as.integer(tau), RLbnd = as.integer(RLbnd), UCL = as.double(UCL), J2 = J2, m = m,
            RL1 = integer(nrep), tauhat = integer(nrep), ocPXL = double(n1 * n2), dQs = double(nrep))
  ocPXL <- matrix(out$ocPXL, ncol = n2, byrow = TRUE)
  return(list(RL1 = out$RL1, tauhat = out$tauhat, ocPXL = ocPXL, dQs = out$dQs))
}


set.seed(314)
system.time(
  GLRrep.out <- GLR_rep(nrep = 1, n1 = n, n2 = n, faultshape = "square", faultsize = fs+0*sqrt(s2), noiselevel = sqrt(s2), tau = tau, RLbnd = 20, 
                        UCL = 12.75, J2 = 4, m = 10)
)

mean(GLRrep.out$RL1)
mean(GLRrep.out$tauhat)
image(GLRrep.out$ocPXL, col = gray(c(0:255)/255))
max(abs(GLRrep.out$ocPXL - PXLoc1st))

betahat <- wav1ds(Z = ims[, , 18], decom.level = J)$BETA
Zhat <- wav1ds(Z = ims[, , 18], decom.level = J)$Zhat
image(ims[, , 18], col = gray(c(0:255)/255))
image(Zhat, col = gray(c(0:255)/255))
img2hat <- wav1ds_img(BETA = betahat, n2 = ncol(ims[, , 2]))
image(img2hat, col = gray(c(0:255)/255))

Dhat <- matrix(0, nrow = n, ncol = n)
Dhat[which(abs(GLRrep.out$ocPXL) > 10^{-15})] <- 1
D <- matrix(0, nrow = n, ncol = n)
D[which(abs(square1) > 10^{-15})] <- 1
image(D, col = gray(c(0:255)/255))
image(Dhat, col = gray(c(0:255)/255))
(0.5 * sum((Dhat - D) > 0)/(n^2 - sum(D)) + (1 - 0.5) * sum((D - Dhat) > 0) / sum(D))
mean(GLRrep.out$dQs)

### Summarize repeated simulations
truetau <- 14
trueRL <- 15

RSalo.CrsFS2 <- read.csv(file = "RSalo_CrsFS2.csv", header = TRUE)
lapply(RSalo.CrsFS2, FUN = mean)
unlist(lapply(RSalo.CrsFS2, FUN = sd))/sqrt(nrow(RSalo.CrsFS2))
lapply(RSalo.CrsFS2, FUN = median)
mean(RSalo.CrsFS2[, "RLs"] == trueRL)
mean(RSalo.CrsFS2[, "taus"] == truetau)

RSalo.CrsFS4 <- read.csv(file = "RSalo_CrsFS4.csv", header = TRUE)
lapply(RSalo.CrsFS4, FUN = mean)
unlist(lapply(RSalo.CrsFS4, FUN = sd))/sqrt(nrow(RSalo.CrsFS4))
lapply(RSalo.CrsFS4, FUN = median)
mean(RSalo.CrsFS4[, "RLs"] == trueRL)
mean(RSalo.CrsFS4[, "taus"] == truetau)

RSalo.SqrFS2 <- read.csv(file = "RSalo_SqrFS2.csv", header = TRUE)
lapply(RSalo.SqrFS2, FUN = mean)
unlist(lapply(RSalo.SqrFS2, FUN = sd))/sqrt(nrow(RSalo.SqrFS2))
lapply(RSalo.SqrFS2, FUN = median)
mean(abs(RSalo.SqrFS2[, "taus"] - 14) < 1/1000)
mean(RSalo.SqrFS2[, "RLs"] == trueRL)
mean(RSalo.SqrFS2[, "taus"] == truetau)


RSalo.SqrFS4 <- read.csv(file = "RSalo_SqrFS4.csv", header = TRUE)
lapply(RSalo.SqrFS4, FUN = mean)
unlist(lapply(RSalo.SqrFS4, FUN = sd))/sqrt(nrow(RSalo.SqrFS4))
lapply(RSalo.SqrFS4, FUN = median)
mean(abs(RSalo.SqrFS4[, "taus"] - 14) < 1/1000)
mean(RSalo.SqrFS4[, "RLs"] == trueRL)
mean(RSalo.SqrFS4[, "taus"] == truetau)


RSscad.CrsFS2 <- read.csv(file = "RSscad_CrsFS2.csv", header = TRUE)
lapply(RSscad.CrsFS2, FUN = mean)
unlist(lapply(RSscad.CrsFS2, FUN = sd))/sqrt(nrow(RSscad.CrsFS2))
lapply(RSscad.CrsFS2, FUN = median)
mean(RSscad.CrsFS2[, "RLs"] == trueRL)
mean(RSscad.CrsFS2[, "taus"] == truetau)

RSscad.CrsFS4 <- read.csv(file = "RSscad_CrsFS4.csv", header = TRUE)
lapply(RSscad.CrsFS4, FUN = mean)
unlist(lapply(RSscad.CrsFS4, FUN = sd))/sqrt(nrow(RSscad.CrsFS4))
lapply(RSscad.CrsFS4, FUN = median)
mean(RSscad.CrsFS4[, "RLs"] == trueRL)
mean(RSscad.CrsFS4[, "taus"] == truetau)

RSscad.SqrFS2 <- read.csv(file = "RSscad_SqrFS2.csv", header = TRUE)
lapply(RSscad.SqrFS2, FUN = mean)
unlist(lapply(RSscad.SqrFS2, FUN = sd))/sqrt(nrow(RSscad.SqrFS2))
lapply(RSscad.SqrFS2, FUN = median)
mean(RSscad.SqrFS2[, "RLs"] == trueRL)
mean(RSscad.SqrFS2[, "taus"] == truetau)

RSscad.SqrFS4 <- read.csv(file = "RSscad_SqrFS4.csv", header = TRUE)
lapply(RSscad.SqrFS4, FUN = mean)
unlist(lapply(RSscad.SqrFS4, FUN = sd))/sqrt(nrow(RSscad.SqrFS4))
lapply(RSscad.SqrFS4, FUN = median)
mean(RSscad.SqrFS4[, "RLs"] == trueRL)
mean(RSscad.SqrFS4[, "taus"] == truetau)


RSglr.CrsFS2 <- read.csv(file = "RSglr_CrsFS2.csv", header = TRUE)
lapply(RSglr.CrsFS2, FUN = mean)
unlist(lapply(RSglr.CrsFS2, FUN = sd))/sqrt(nrow(RSglr.CrsFS2))
lapply(RSglr.CrsFS2, FUN = median)

RSglr.CrsFS4 <- read.csv(file = "RSglr_CrsFS4.csv", header = TRUE)
lapply(RSglr.CrsFS4, FUN = mean)
unlist(lapply(RSglr.CrsFS4, FUN = sd))/sqrt(nrow(RSglr.CrsFS4))
lapply(RSglr.CrsFS4, FUN = median)

RSglr.SqrFS2 <- read.csv(file = "RSglr_SqrFS2.csv", header = TRUE)
lapply(RSglr.SqrFS2, FUN = mean)
unlist(lapply(RSglr.SqrFS2, FUN = sd))/sqrt(nrow(RSglr.SqrFS2))
lapply(RSglr.SqrFS2, FUN = median)
mean(abs(RSglr.SqrFS2[, "taus"] - 14) < 1/1000)


RSglr.SqrFS4 <- read.csv(file = "RSglr_SqrFS4.csv", header = TRUE)
lapply(RSglr.SqrFS4, FUN = mean)
unlist(lapply(RSglr.SqrFS4, FUN = sd))/sqrt(nrow(RSglr.SqrFS4))
lapply(RSglr.SqrFS4, FUN = median)
mean(abs(RSglr.SqrFS4[, "taus"] - 14) < 1/1000)


## Tile images

library(jpeg)
tile <- readJPEG(source = "Tile_Nom.jpeg")
tile <- t(tile)[, ncol(tile):1]
image(tile, col = gray.colors(256))

n1 <- nrow(tile)
n2 <- ncol(tile)
tile_fault <- matrix(0, nrow = n1, ncol = n2)
for (k1 in 1:n1) {
  for (k2 in 1:n2) {
    x1 <- ((k1 - 1) / n1 - 0.85)
    x2 <- ((k2 - 1) / n2 - 0.6)
    xx1 <- (x1 - x2) / sqrt(2)
    xx2 <- (x1 + x2) / sqrt(2)
    if ((xx1/0.05)^2 + (xx2/0.005)^2 <= 1) {
      tile_fault[k1, k2] <- 0.015
    }
  }
}
par(mfrow = c(1, 1))
image(tile + tile_fault, col = gray.colors(256)) 

s.tile <- 0.01
n_ims.tile <- 30
tau.tile <- 20

ims.tile <- array(NA, dim = c(n1, n2, n_ims.tile))
set.seed(314)
for (t in 1:n_ims.tile) {
  if (t < tau.tile) {
    ims.tile[, , t] <- tile + matrix(rnorm(n1 * n2, sd = s.tile), ncol = n2, byrow = TRUE)  
  } else {
    ims.tile[, , t] <- tile + tile_fault + matrix(rnorm(n1 * n2, sd = s.tile), ncol = n2, byrow = TRUE)  
  }
}

par(mfrow = c(2, 2))
for (t in (tau.tile - 1):(tau.tile + 1)) {
  image(ims.tile[, , t], col = gray.colors(256), zlim = range(ims.tile))
}

### ALASSO: tile image 

beta_tile <- wav2d_coef(Z = tile)$wcoef
p.tile <- length(beta_tile)

### IC stats for tile image: run by ICstats_Tile.R

nrep <- 25000
ncores <- 10
nrpc <- as.integer(nrep / ncores)
RNGkind("L'Ecuyer-CMRG")
set.seed(314)
system.time(
  alo.mcout <- mcmapply(ALO_ICstat, n_IC = rep(nrpc, ncores), MoreArgs = list(p = p.tile, sig2 = s.tile^2/(n1 * n2)), mc.cores = ncores)
)
invisible(rnorm(100)) # reset the random seed

alo.bar_mat <- alo.mcout[1:p, ]
alo.S2_mat <- alo.mcout[(p+1): (2*p), ]^2
alo.ic.bar <- rowMeans(alo.bar_mat)
alo.ic.std <- sqrt(rowMeans(alo.S2_mat) + rowMeans(alo.bar_mat^2) - alo.ic.bar^2)

### Read in output from ICstats_Tile.R

ICstats_Tile <- read.csv(file = "ICstats_Tile.csv")

alo.ic.bar.tile <- ICstats_Tile[, "alo.ic.bar"]
alo.ic.std.tile <- ICstats_Tile[, "alo.ic.std"]
scad.ic.bar.tile <- ICstats_Tile[, "scad.ic.bar"]
scad.ic.std.tile <- ICstats_Tile[, "scad.ic.std"]

### ALASSO Tile control limit: 4.065

ucl.alo <- 4.065
alo.Qs.tile <- ALO_ims(imseq = ims.tile, wcoef_mean = beta_tile, sig2 = s.tile^2/(n1 * n2), IC.bar = alo.ic.bar.tile, IC.std = alo.ic.std.tile, 
                  la = 0.2, q = p.tile)
(OCindx <- which(alo.Qs.tile >= ucl.alo, arr.ind = TRUE))
(tauhat <- change_loc(ims = ims.tile[, , 1:OCindx[1]], wcoef_mean = beta_tile))
alo.OC.out <- ALO_change_wcoef(ims.OC = ims.tile[, , (tauhat$loc + 1):OCindx[1]], sig2 = s.tile^2/(n1 * n2), wcoef_mean = beta_tile)
alo.img.diag <- wav2d_img(n1 = n1, n2 = n2, J1 = as.integer(log(n1, base = 2)), J2 = as.integer(log(n2, base = 2)), wcoef = alo.OC.out$OC_wcoef)$img
image(tile + alo.img.diag, col = gray.colors(256))
image(tile + tile_fault, col = gray.colors(256))
which.min(alo.OC.out$EBIC)
which(abs(alo.OC.out$OC_wcoef) > 0)
par(mfrow = c(2, 2))
image(abs(alo.img.diag) > 0, col = gray.colors(256))
image(tile_fault, col = gray.colors(256))
par(mfrow = c(1, 1))
image(tile + 20 * tile_fault + (abs(alo.img.diag) > 0), col = gray.colors(256))
image(tile + 20 * tile_fault, col = gray.colors(256))
rect(xleft = 0.81, ybottom = 0.59, xright = 0.88, ytop = 0.63, border = "red", lty = "dashed", lwd = 2)

plot(alo.Qs.tile, type = "b", ylim = range(c(alo.Qs.tile, ucl.alo + 1)))
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.alo, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")


### SCAD Tile control limit: 4.037

ucl.scad <- 4.021
scad.Qs.tile <- SCAD_ims(imseq = ims.tile, wcoef_mean = beta_tile, sig2 = s.tile^2/(n1 * n2), IC.bar = scad.ic.bar.tile, IC.std = scad.ic.std.tile, 
                       la = 0.2, q = p.tile)
(OCindx.scad <- which(scad.Qs.tile >= ucl.scad, arr.ind = TRUE))
(tauhat.scad <- change_loc(ims = ims.tile[, , 1:OCindx.scad[1]], wcoef_mean = beta_tile))
scad.OC.out <- SCAD_change_wcoef(ims.OC = ims.tile[, , (tauhat.scad$loc + 1):OCindx.scad[1]], sig2 = s.tile^2/(n1 * n2), wcoef_mean = beta_tile)
scad.img.diag <- wav2d_img(n1 = n1, n2 = n2, J1 = as.integer(log(n1, base = 2)), J2 = as.integer(log(n2, base = 2)), wcoef = scad.OC.out$OC_wcoef)$img
image(tile + scad.img.diag, col = gray.colors(256))
image(tile + tile_fault, col = gray.colors(256))
which.min(scad.OC.out$EBIC)
which(abs(scad.OC.out$OC_wcoef) > 0)
par(mfrow = c(2, 2))
image(abs(scad.img.diag) > 0, col = gray.colors(256))
image(tile_fault, col = gray.colors(256))
par(mfrow = c(1, 1))
image(tile + 20 * tile_fault + (abs(scad.img.diag) > 0), col = gray.colors(256))
image(tile + 20 * tile_fault, col = gray.colors(256))
rect(xleft = 0.81, ybottom = 0.59, xright = 0.88, ytop = 0.63, border = "red", lty = "dashed", lwd = 2)

plot(scad.Qs.tile, type = "b", ylim = range(c(scad.Qs.tile, ucl.scad + 1)))
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.scad, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")


### GLR Tile control limit

m <- 10
nrpercore <- 500
ncore <- 10
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(314) # takes 5 minutes to run
system.time(
  IdpRLs <- mcmapply(IdpGLR_RLs0, nRuns = rep(nrpercore, ncore), MoreArgs = list(nGLRs = n2 * 2^4, upbnd = 200 * 10, wsize = m,
                                                                                 UCL = 13.444),
                     mc.cores = detectCores())
)
invisible(rnorm(100)) # reset seed
mean(IdpRLs)
sd(IdpRLs)/sqrt(nrpercore * ncore)

ucl.glr <- 13.444
J <- 4
betastack_mean <- c(t(wav1ds(Z = tile, decom.level = J)$BETA))
GLRs.out <- GLR_ims(imseq = ims.tile, decom.level = J, betasig2 = s.tile^2/n1, m = 10, wcoef_mean = betastack_mean)

plot(GLRs.out$GLRchart, type = "b", ylim = range(c(GLRs.out$GLRchart, ucl.glr + 1)))
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.glr, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")

### Make figure for publication
layout(matrix(c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 4), rep(6, 4), rep(7, 4)), 2, 12, byrow = TRUE))
par(mar = c(4, 3, 2, 3))
image(ims.tile[, , 1], col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(a)", line = 0.5)
image(ims.tile[, , 30] + 19 * tile_fault, col = gray.colors(256), xaxt = "n", yaxt = "n")
title(sub = "(b)", line = 0.5)
image(tile + 20 * tile_fault, col = gray.colors(256), xaxt = "n", yaxt = "n")
rect(xleft = 0.81, ybottom = 0.59, xright = 0.88, ytop = 0.63, border = "red", lty = "solid", lwd = 2)
title(sub = "(c)", line = 0.5)
image(tile + 20 * tile_fault, col = gray.colors(256), xaxt = "n", yaxt = "n")
rect(xleft = 0.81, ybottom = 0.59, xright = 0.88, ytop = 0.63, border = "red", lty = "solid", lwd = 2)
title(sub = "(d)", line = 0.5)

plot(alo.Qs.tile, type = "b", ylim = range(c(alo.Qs.tile, ucl.alo + 1)), axes = FALSE, xlab = "", ylab = "")
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.alo, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4)
title(xlab = "Time Index", line = 1.5)
title( ylab = expression(Q["AL, t"]), line = 1.5, cex = 1.9)
title(sub = "(e)", line = 2.3)

plot(scad.Qs.tile, type = "b", ylim = range(c(scad.Qs.tile, ucl.scad + 1)), axes = FALSE, xlab = "", ylab = "")
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.scad, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4)
title(xlab = "Time Index", line = 1.5)
title( ylab = expression(Q["SCAD, t"]), line = 1.5, cex = 1.9)
title(sub = "(f)", line = 2.3)

plot(GLRs.out$GLRchart, type = "b", ylim = range(c(GLRs.out$GLRchart, ucl.glr + 1)), axes = FALSE, xlab = "", ylab = "")
lines(x = 0:(n_ims.tile + 1), y = rep(ucl.glr, n_ims.tile + 2), lty = "dashed", lwd = 2, col = "blue")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4)
title(xlab = "Time Index", line = 1.5)
title( ylab = "KNM Charting Stat.", line = 1.5, cex = 1.2)
title(sub = "(g)", line = 2.3)
