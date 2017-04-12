# Distribuição Katz Lagrangiana.

Versão totalmente draft das funções d, p, q e r para a LKD em R puro.

## Instalação do pacote
```{R}
devtools::install_github('evandeilton/lkd')
```
## Testes
```{R}

# Carrega pacotes para testes

require(data.table)
require(bbmle)
require(fitdistrplus)

# Simula dados da LKD
dados <- rLKD(300, 5, 0.05, 0.1)

# A LKD é muito sensível ao seu espaço paramtértico. Este grid é uma tentativa de restringir o espaço de busca para os estimadores de Máxima Verossimilhança da LKD.

grid <- as.data.table(expand.grid(a = seq(0.01, 10, l = 20), b = seq(0.002, 10, l = 20), beta = seq(-0.001, 1, l = 20)))
grid$valida <- ifelse(grid$a > 0 & grid$b > -1 & grid$beta > -grid$b & grid$beta < 1, 1, 0)
grid <- subset(grid, grid$valida == 1, select = -4)
grid <- grid[sample(nrow(grid))]

# Limites inferior, superior e iniciais para "L-BFGS-B"
lo <- round(sapply(grid, min), 4)
up <- round(sapply(grid, max), 4)
ii <- sapply(grid, mean)

# Looping com fitdist busancoa mle's
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

# Vizualização do ajuste, se houver.
plot(fit)
cdfcomp(fit, addlegend=FALSE)
denscomp(fit, addlegend=FALSE, xlim = c(0,1))
ppcomp(fit, addlegend=FALSE, xlim = c(0,1))

# Análise com mle2 (pacote bbmle)
ll2 <- function(y, a, b, beta){
  l <- NA
  
  if(a > 0 & b > -1 & beta > -b & beta < 1) {
    l <- dLKD(x = y, a=a, b=b, beta=beta, log=TRUE)  
    #j <- !is.finite(l)
    #l[j] <- log(.Machine$double.xmin * (1.15e-16))
    l <- -sum(l, na.rm = T)
  }
  return(l)
}

fit <- mle2(ll2,
            start = as.list(ii), 
            method = "L-BFGS-B",
            optimizer = "optim",
            lower = as.list(lo),
            upper = as.list(up),
            data = list(y = dados),
            control = list(fnscale = -1, trace = T, maxit = 10000, factr = 1e-10),
            use.ginv = TRUE)
summary(fit)
pro <- profile(fit)
plot(pro)
```
