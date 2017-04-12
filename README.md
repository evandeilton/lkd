# Distribuição Katz Lagrangiana.

Versão totalmente draft das funções d, p, q e r para a LKD em R puro.

## Instalação do pacote
```{R}
devtools::install_github('evandeilton/lkd')
```
## Testes
```{R}
################## Análises por MLE ######################

## Pacotes
require(fitdistrplus)
require(lkd)
require(bbmle)
require(data.table)

## Dados simulados
set.seed(2)
N <- 2000
x <- 0:(N-1)
a <- 3; b<- 0.5; beta <- 0.05
dados <- rLKD(N, a, b, beta)

par(mfrow = c(2,2))
hist(dados, nclass = 8, prob=T,main="Amostra",cex.axis=1,cex.lab=1,col="darkgreen")
curve(dLKD(x, a, b, beta), 0, max(y)*2, typ='h', main="Distribuição",cex.axis=1,cex.lab=1,col="darkgreen")
plot(function(x) pLKD(x, a, b, beta), 0, max(y)*2, main = "Cumulativa")
xx <- pLKD(0:80, a, b, beta)
plot(function(x) qLKD(x, a, b, beta), min(xx), max(xx), main = "Inversa Cumulativa")

## Grid para gerar intervalo de busca. A LKD é sensível ao seu espaço paramétrico.
grid <- as.data.table(expand.grid(a = seq(0.01, 10, l = 20), b = seq(0.002, 10, l = 20), beta = seq(-0.2, 1, l = 20)))
grid$valida <- ifelse(grid$a > 0 & grid$b > -1 & grid$beta > -grid$b & grid$beta < 1, 1, 0)
grid <- subset(grid, grid$valida == 1, select = -4)
grid <- grid[sample(nrow(grid))]

## Valores para limites inferior, superior e chutes iniciais para o L-BFGS-B
lo <- round(sapply(grid, min), 4)
up <- round(sapply(grid, max), 4)
ii <- sapply(grid, mean)

## Looping usando a função fitdist
i <- 1
while(i <= 10){
  ini <- c(a = as.numeric(grid[i, 1]), b=as.numeric(grid[i, 2]), beta = as.numeric(grid[i, 3]))
  fit <- try(fitdist(dados, "LKD", method = "mle", start = ini, lower = lo, upper = up, discrete = TRUE))
  if(class(fit)[1] != "try-error"){
    return(fit) 
    break      
  } else {
    cat("Erro em:", unlist(ini), "\n") 
  }
  i <- i + 1
}

## Visualização do ajuste dos MLE's
summary(fit)
plot(fit)
cdfcomp(fit, addlegend=FALSE)
denscomp(fit, addlegend=FALSE, xlim = c(0,1))
ppcomp(fit, addlegend=FALSE, xlim = c(0,1))


## Log-verossimilhança usando a densidade LKD para a função mle2
ll2 <- function(y, a, b, beta){
  l <- NA
  
  if(a > 0 & b > -1 & beta > -b & beta < 1) {
    l <- dLKD(x = y, a, b, beta, log=TRUE)  
    #j <- !is.finite(l)
    #l[j] <- log(.Machine$double.xmin * (1.15e-16))
    l <- -sum(l, na.rm = T)
  }
  return(l)
}


## Log-verossimilhança usando a definição de Consul e Famoye em no capítulo Lagrangia Katz Distribution (Prem C. Consul, Felix Famoye, Samuel Kotz-Lagrangian Probability Distributions-Birkhäuser (2006))


ll <- function(par, y){
  #a, b, beta
  a <- par[1];b<-par[2]; beta <- par[3]
  xbar <- mean(y)
  sbar <- sd(y)
  n  <- length(y)
  n0 <- sum(y == y[1])
  
  ni <- table(y)
  s1 <- sapply(2:length(ni), function(i) {
    ni[i]*log(i)
  })
  
  #ii <- !is.finite(s1)
  #s1[ii] <- log(.Machine$double.xmin * (1.15e-16))
  s1 <- sum(unlist(s1), na.rm = T)
  
  s2 <- sapply(2:length(ni), function(i) {
    sapply(1:(i-1), function(j) {
      ni[i]*log(a + b*i + beta*j)
    })
  })
  s2 <- unlist(s2)
  
  #jj <- !is.finite(s2)
  #s2[jj] <- log(.Machine$double.xmin * (1.15e-16))
  s2 <- sum(s2, na.rm = T)
  
  #cat("Soma:", -su, "\n")
  
  l <- (n-n0)*log(a) + (n/beta)*(a + b*xbar)*log(1-beta) - s1 + s2
  
  return(-sum(l, na.rm = T))
}


## Ajuste por bbmle. Note que fnscale = 1, pois a função de verossimilhança já vem negativa.
parnames(ll) <- c("a", "b", "beta")
fit2 <- mle2(ll,
             vecpar = TRUE,
             method = "L-BFGS-B",
             optimizer = "optim",             
             start = as.list(ii), 
             lower = as.list(lo),
             upper = as.list(up),
             data = list(y = dados),
             control = list(fnscale = 1, trace = T, maxit = 10000, factr = 1e-10),
             use.ginv = TRUE
)
summary(fit2)
coef(fit2)
vcov(fit2)
pr <- profile(fit2)
plot(pr)


y <- dados
LL <- Vectorize(
  FUN=function(va, vb, vbeta, y){
    ll(c(va, vb, vbeta), y=y)
  },
  vectorize.args=c("va", "vb","vbeta"))

## Grid para gerar intervalo de busca. A LKD é sensível ao seu espaço paramétrico.
grid <- as.data.table(expand.grid(a = seq(0.01, 10, l = 30), b = seq(0.002, 10, l = 30), beta = seq(-0.2, 1, l = 20)))
grid$valida <- ifelse(grid$a > 0 & grid$b > -1 & grid$beta > -grid$b & grid$beta < 1, 1, 0)
grid <- subset(grid, grid$valida == 1, select = -4)

pars <- coef(fit2)
ci <- confint(fit2, method="quad")
logL <- as.numeric(logLik(fit2))

va <- seq(min(grid$a), 5, l=100)
vb <- seq(min(grid$b), 2, l=100)
vbeta <- seq(-1, max(grid$beta), l=100)

lla <- outer(va, vb, vbeta, FUN = LL, y = y)

fnPlot <- function(y, vx, vy, cx, cy, cix, ciy, npar, logL){
  contour(vx, vy, y,
          xlab = names(cx),
          ylab = names(cy))
  
  abline(v=cx, h=cy, lty=2)
  
  points(cx, cy, type = "o", col = "red", cex = 1.5)
  text(cx, cy*1.2, labels = paste("logL:", round(logL, 2)))
  
  contour(vx, vy, y, add=TRUE,
          levels=(logL-0.5*qchisq(0.95, df=npar)),
          col=2)
  
  abline(v=cix, h=ciy, lty=3, col=2)
}

par(mfrow = c(1,3))
fnPlot(lla, va, vb, pars[1], pars[2], ci[1], ci[2], 2, logL)
fnPlot(lla, va, vbeta, pars[1], pars[3], ci[1], ci[3], 2, logL)
fnPlot(lla, vb, vbeta, pars[2], pars[3], ci[2], ci[3], 2, logL)
```
