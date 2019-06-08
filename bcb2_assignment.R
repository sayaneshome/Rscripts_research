Q <- matrix(rgamma(16, shape=1), nrow=4)
diag(Q) <- 0
diag(Q) <- -rowSums(Q)
library(Matrix)
# compute P(1)
expm(Q)
# estimate stationary distribution
pi <- expm(Q*10)[1,]
# scale Q so that time is measured as "evolutionary time"
Q <- Q/sum(- pi * diag(Q))
# compute P(1)
expm(Q)

1 - sum(pi * diag(expm(Q)))
## [1] 0.4378962
# PAM1 is P(0.01)
PAM1 <- expm(Q*0.01)
# in discrete time, step forward 32 steps, or approximately
# 0.32 units of evolutionary time
PAM2 <- PAM1 %*% PAM1
PAM4 <- PAM2 %*% PAM2
PAM8 <- PAM4 %*% PAM4
PAM16 <- PAM8 %*% PAM8
PAM32 <- PAM16 %*% PAM16
PAM32

S <- round(log(t(t(PAM32)/pi)))
S
## 4 x 4 Matrix of class "dgeMatrix"
## [,1] [,2] [,3] [,4]
## [1,] 2 -1 1 -2
## [2,] 0 1 -2 -1
## [3,] 0 -1 2 -1
## [4,] -1 -1 -1 0
sum(diag(pi) %*% S %*% diag(pi))


r <- c(3.821, 21.193, 2.615, 2.567, 34.546, 1.000)
pi <- c(0.298, 0.249, 0.232, 0.221)
R <- matrix(c(0, r[c(1:3, 1)], 0, r[c(4:5,2,4)], 0, r[c(6,3,5,6)], 0),byrow=T, ncol=4)  # change here
Q <- R %*% diag(pi)
diag(Q) <- -rowSums(Q)
P <- diag(1/rowSums(Q - diag(diag(Q)))) %*% Q
set.seed(1111)
nucs <- c('A','C','G','T')
#simulate 50,100 and 1000 sequence alignment

Q.norm <- -Q/sum(pi * diag(Q))
P.norm <- diag(1/rowSums(Q.norm - diag(diag(Q.norm)))) %*% Q.norm
v.norm <- -diag(Q.norm)
diag(P.norm) <- rep(0, 4)
t <- 0
t.changes <- x <- NULL
x[1] <- sample(1:4, size=1, prob=pi)
i <- 1
while (t < 2.5)
{#t.changes[i] <- 2.5*(t + rexp(1, rate=v.norm[x[i]]))
  t.changes[i] <- t + rexp(1, rate=v.norm[x[i]])
  t <- t.changes[i]
  if (t >= 2.5) break
  x[i+1] <- sample(1:4, size=1, prob=P.norm[x[i],])
  i <- i + 1
}

library(expm, quietly = T, warn.conflicts = F)
P.1 <- expm(Q.norm*0.1)
print(P.1)
P.1 <- P.1 %^% 10
print(P.1)
P.1 <- expm(Q.norm*0.01)
P.10 <- P.1 %^% 10

S <- log(P.10 %*% diag(1/pi))
for (k in c(1,20*(1:10)))
  {S.tmp <- round(k * S)
  err <- sum(abs(S.tmp - k * S))/sum(k*S)
  s <- sum(diag(pi) %*% S.tmp %*% diag(pi))
  cat(k, s, err, "\n")
  # multiplier, avg. score, rel. error
}
k.choose <- 100
S.tmp <- round(k.choose * S)

f <- function(t, pi, S)
{sum(diag(pi) %*% exp(t*S) %*% diag(pi)) - 1}
lambda <- uniroot(f, interval=c(1e-5,2), pi=pi, S=S.tmp)$root
S <- as.vector(S.tmp)# the scoring matrix
prob <- as.vector(pi %*% t(pi))# score probabilities under H0
n.rep <- 100000# no. of simulation replicates
H <- rep(0, n.rep)# excursion heights
L <- rep(-1, n.rep)# excursion lengths
Qt <- rep(0, max(S))# prob first positive integer visited
R <- rep(0, -min(S))# prob score at excursion 

score <- 0# current excursion score, initial 0
i <- 1
while (i < n.rep){
  # sample score under H0
  score <- score + sample(S, size=1, prob=prob)
  L[i] <- L[i] + 1# and length# first positive value visited
  if (score > 0 & L[i] == 0)
  {Qt[score] <- Qt[score] + 1}# record maximal height
  if (H[i] < score) H[i] <- score# if score drops below 0, the excursion ends: reset
  if (score < 0){
    # record score at excursion end
    R[-score] <- R[-score] + 1
    i <- i + 1# next excursion
    score <- 0 # reset score to 0
  }}

R <- R/n.rep
Qt <- Qt/n.rep
Q.bar <- 1 - sum(Qt, na.rm=T)# prob. no positive value visited
C <- Q.bar*(1 - sum(R * exp(-lambda * 1:length(R)), na.rm=T)) /(1 - exp(-lambda)) /sum(1:length(Qt)*Qt*exp(1:length(Qt)*lambda), na.rm=T)
A <- -sum(1:length(R)*R, na.rm=T)/sum(S * prob)
