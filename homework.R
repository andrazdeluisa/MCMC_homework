library(mcmcse)
library(ggplot2)
library(mvtnorm)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### EXERCISE 1

set.seed(0)

f <- function(x) 10 * exp(-5 * (x - 3) ** 4)
p <- function(x) dnorm(x, mean=0, sd=1, log=FALSE)
g <- function(x) f(x) * p(x)

# quadrature method
I1 <- integrate(g, -Inf, Inf)$value

# monte carlo
x <- rnorm(100, mean=0, sd=1)
fx <- f(x)
V2 <- var(fx)
I2 <- mean(fx)

# repeated monte carlo
I3 <- numeric()
ci <- numeric()
for (i in 1:1000){
  x_tmp <- rnorm(100, mean=0, sd=1)
  fx_tmp <- f(x_tmp)
  I3[i] <- mean(fx_tmp)
  s <- sd(fx_tmp)
  err <- s * qnorm(0.95)/sqrt(100)
  ci[i] <- (I1 > I3[i] - err) & (I1 < I3[i] + err)
}

mean(ci)

mean(I3)
err <- sd(I3) * qnorm(0.95)/sqrt(1000)
mean(I3) - err
I1
mean(I3) + err

# importance sampling
x2 <- rnorm(100, mean=3, sd=0.8)
h <- function(x) g(x) / dnorm(x, mean=2.5, sd=.6)
f2x <- h(x2)
I4 <- mean(f2x)
V4 <- var(f2x)


xx <- seq(1, 4, 0.005)
plot(xx, h(xx), col='red', type='l')
points(xx, g(xx), type='l')

# rejection sampling

M <- 4 / sqrt(2*pi)
i <- 1
env <- function(x) dlogis(x, location=0, scale=1)
x3 <- numeric()

while (i <= 100){
  xi <- rlogis(1, location=0, scale=1)
  ui <- runif(1, 0, M * env(xi))
  if (ui < p(xi)){
    x3[i] <- xi
    i <- i + 1
  }
}

fx3 <- f(x3)
V5 <- var(fx3)
I5 <- mean(fx3)

# metropolis hastings
chain_rule1 <- function(x, delta){
  runif(1, x - delta, x + delta)
}

mh <- function(x0, n, proposal, p){
  x <- matrix(nrow=n, ncol=length(x0))
  xi <- x0
  rej <- 0
  for (i in 1:n){
    xj <- proposal(xi)
    alpha <- min(1, p(xj) / p(xi))
    if (p(xi) == 0) alpha <- 1
    u <- runif(1, 0, 1)
    if (u < alpha){
      xi <- xj
      x[i,] <- xj
    } else {
      x[i,] <- xi
      rej <- rej + 1
    }
  }
  return(list(x=x, rejection_rate=rej/n))
}

x4 <- mh(0, 100, function(x) chain_rule1(x, 2), p)
x4$rejection_rate
fx4 <- f(as.numeric(x4$x))

I6 <- mean(fx4)
V6 <- var(fx4)
ess(x4$x)
acf(x4$x)

### EXERCISE 3

set.seed(1)

f_banana <- function(x) {
  exp(-(x[1]^2)/200- 0.5 * (x[2]+ 0.05 * x[1]^2 - 100*0.05)^2 )
}

chain_rule2 <- function(x){
  y1 <- runif(1, x[1] - 1, x[1] + 1)
  y2 <- runif(1, x[2] - 1, x[2] + 1)
  y <- c(y1, y2)
  return(y)
}

chain_rule3 <- function(x){
  y1 <- runif(1, x[1] - 20, x[1] + 20)
  y2 <- runif(1, x[2] - 20, x[2] + 20)
  y <- c(y1, y2)
  return(y)
}

chain_rule4 <- function(x){
  y <- rmvnorm(1, mean=x, sigma=diag(c(16, 25)))
  return(y)
}

chain_rule_best <- function(x){
  y1 <- runif(1, x[1] - 5, x[1] + 5)
  y2 <- runif(1, x[2] - 5, x[2] + 5)
  y <- c(y1, y2)
  return(y)
}

xs <- seq(-25, 25, 0.5)
ys <- seq(-25, 25, 0.5)

z <- matrix(nrow=length(xs), ncol=length(ys))
xxs <- matrix(nrow=length(xs), ncol=length(ys))
yys <- matrix(nrow=length(xs), ncol=length(ys))
for (i in 1:length(xs)){
  for (j in 1:length(ys)){
    z[i, j] <- f_banana(c(xs[i], ys[j]))
    xxs[i, j] <- xs[i]
    yys[i, j] <- ys[j]
  }
}

df <- data.frame(x=as.numeric(xxs), y=as.numeric(yys), z=as.numeric(z))

starting_points <- matrix(c(0, 20, -10, 5, -15, -10), nrow=3, ncol=2)
rules <- list(chain_rule2, chain_rule3, chain_rule_best)

chains <- list()
rejections <- list()
essizes <- list()

for (i in 1:3){
  chains_tmp <- list()
  rej_tmp <- c()
  ess_tmp <- matrix(nrow=3, ncol=2)
  for (j in 1:3){
    tmp <- mh(starting_points[j,], 1000, rules[[i]], f_banana)
    chains_tmp[[j]] <- tmp$x
    rej_tmp[j] <- tmp$rejection_rate
    ess_tmp[j,] <- ess(tmp$x)
  }
  chains[[i]] <- chains_tmp
  rejections[[i]] <- rej_tmp
  essizes[[i]] <- ess_tmp
  
  ggplot(df, aes(x=x, y=y)) + geom_contour(aes(z=z)) + 
          geom_jitter(data=as.data.frame(chains_tmp[[1]][1:100,]), aes(x=V1, y=V2, z=0), col='red', size=.1) +
          geom_line(data=as.data.frame(chains_tmp[[1]][1:100,]), aes(x=V1, y=V2, z=0), col='red', size=.05) +
          geom_line(data=as.data.frame(chains_tmp[[2]][1:100,]), aes(x=V1, y=V2, z=0), col='green', size=.05) +
          geom_jitter(data=as.data.frame(chains_tmp[[2]][1:100,]), aes(x=V1, y=V2, z=0), col='green', size=.1) +
          geom_line(data=as.data.frame(chains_tmp[[3]][1:100,]), aes(x=V1, y=V2, z=0), col='black', size=.05) +
          geom_jitter(data=as.data.frame(chains_tmp[[3]][1:100,]), aes(x=V1, y=V2, z=0), col='black', size=.1) +
          theme_bw() +
          ggsave(paste0('rule', i, '.pdf', collapse=''), dpi=600, width=5, height=5)
}

rejections
essizes

### EXERCISE 4

set.seed(1)
x <- rnorm(500, mean=0, sd=1)
x <- matrix(x, nrow=100, ncol=5)
df_toy <- as.data.frame(x)
colnames(df_toy) <- c("x1", "x2", "x3", "x4", "x5")
df_toy$p <- 1 / (1 + exp(-3*df_toy$x1 - 2*df_toy$x2 + df_toy$x3/10 + 1))
df_toy$y <- rbinom(100, 1, df_toy$p)
df_toy$p <- NULL

sigmoid <- function(x){
  1 / (1 + exp(-x))
}

logreg_lik <- function(x, y, beta){
  interc <- beta[length(beta)]
  beta2 <- beta[1:(length(beta) - 1)]
  mu <- sigmoid(x %*% beta2 + interc * rep(1, length(x[,1])))
  tmp <- mu ** y * (1 - mu) ** (1 - y)
  return(prod(tmp))
}

posterior <- function(x, y, beta){
  logreg_lik(x, y, beta) * prod(dnorm(beta, mean=0, sd=10))
}

post_distr <- function(beta) posterior(as.matrix(df_toy[c("x1", "x2", "x3", "x4", "x5")]), as.matrix(df_toy$y), beta)

chain_rule5 <- function(x){
  y <- rmvnorm(1, mean=x, sigma=diag(rep(1, length(x))))
  return(y)
}

starting_points <- matrix(rnorm(24, mean=0, sd=2), nrow=4, ncol=6)

chains_tmp <- list()
rej_tmp <- c()
ess_tmp <- matrix(nrow=4, ncol=6)
for (j in 1:4){
  tmp <- mh(starting_points[j,], 10000, chain_rule5, post_distr)
  chains_tmp[[j]] <- tmp$x
  rej_tmp[j] <- tmp$rejection_rate
  ess_tmp[j,] <- ess(tmp$x)
}

rej_tmp
ess_tmp

acf(chains_tmp[[3]][,1])
acf(chains_tmp[[1]][,2])
acf(chains_tmp[[2]][,3])
acf(chains_tmp[[2]][,4])
acf(chains_tmp[[3]][,5])
acf(chains_tmp[[3]][,6])

df <- as.data.frame(chains_tmp[[3]])
df$ID <- 1:10000

# traceplot (play with y=V1, ... V6 for changing between variables)
g1 <- ggplot(df, aes(x = ID, y = V1)) + 
  geom_line()
plot(g1)

# improved proposal

sigma1 <- matrix(nrow=6, ncol=6)
sigma1[1,] <- c(20, -10, 1, 0, 0, 0)
sigma1[2,] <- c(-10, 15, 0, 0, 0, 0)
sigma1[3,] <- c(1, 0, 2, 0, 0, 0)
sigma1[4,] <- rep(0, 6)
sigma1[5,] <- rep(0, 6)
sigma1[6,] <- c(0, 0, 0, 0, 0, 1)
sigma2 <- diag(c(10, 10, 5, 1, 1, 2))
sigma3 <- diag(c(.5, .5, .5, .1, .1, .5))
sigma4 <- diag(rep(.1, 6))
sigma5 <- diag(rep(.4, 6))

sigma <- sigma3

chain_rule6 <- function(x){
  y <- rmvnorm(1, mean=x, sigma=sigma)
  return(y)
}

tmp <- mh(starting_points[1,], 10000, chain_rule6, post_distr)
chains_tmp2 <- tmp$x
rej_tmp2 <- tmp$rejection_rate
ess_tmp2 <- ess(tmp$x)

rej_tmp2
ess_tmp2
#mcse.mat(chains_tmp[[2]])

acf(chains_tmp2[,1])
acf(chains_tmp2[,2])
acf(chains_tmp2[,3])
acf(chains_tmp2[,4])
acf(chains_tmp2[,5])
acf(chains_tmp2[,6])

df2 <- as.data.frame(chains_tmp2)
df2$ID <- 1:10000

# traceplot (play with y=V1, ... V6 for changing between variables)
g2 <- ggplot(df2, aes(x = ID, y = V1)) + 
  geom_line()
plot(g2)

# posterior means of parameters
beta_est <- mcse.mat(chains_tmp2)
beta_est

# posterior probability
mean(abs(df2$V3) > .1)
