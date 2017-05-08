#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Geração de números aleatórios da LKD pelo método da inversa da função de probabilidade acumulada.
#' @param n número de observacções
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param ... passagem de argumentos
#'
#' @examples rLKD(50, a = 5, b = 0.05, beta = 0.1)
#' @export
rLKD <- function(n, a = 5, b = 0.02, beta = 0.3, ...){
  ve <- sort(runif(n, min = 0, max = 1))
  oo <- c()

  for(i in seq_along(ve)){
    o <- .LKG_one(u = ve[i], a = a, b = b, beta = beta)
    oo[[i]] <- o
  }

  return(sample(oo, replace = FALSE))
}

## @export
.LKG_one <- function(u, a = 5, b = 0.02, beta = 0.3, r = TRUE, ...){
  X <- 0
  W <- (1-beta)^(a/beta)
  S <- P <- W
  U <- u #sort(runif(1, min = 0, max = 1))
  if(U <= S){
    if(!r) return(S)
    else return(X) #cat("x:",X)
  }
  X <- X + 1
  P <- S * a*(1-beta)^(b/beta)
  S <- S + P
  if(U <= S) {
    if(!r) return(S)
    else return(X) #cat("x:",X)
  }
  while(U > S){
    X <- X + 1
    k <- 1:(length(1:X)-1)
    C <- prod(1 + b/(a+b*X+beta*k))
    P <- ((a + b*(X+1) + beta*X)/(X+1))*(1-beta)^(b/beta)*C*S
    S <- S + P

    ##cat("x: ", X, "u: ", U, "c: ", C, "s: ", S, "\n")
  }
  if(!r){
    return(S)
  } else return(X) #cat("x:",X)
}



#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Probabilidades acumuladas da LKD.
#' @param q quantis de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param log se TRUE retorna o log das probabilidades
#' @param ... passagem de argumentos
#' @details todo
#'
#' @examples pLKD(q = 50, a = 5, b = 0.05, beta = 0.1)
#' @export
pLKD <- function(q, a = 5, b = 0.02, beta = 0.3, log = F, ...){
  #Q <- q
  #if(Q == 0) Q <- 1
  #ve <- sort(runif(Q, min = 0, max = 1))
  #oo <- c()

  #for(i in seq_along(ve)){
  #  o <- .LKG_one(u = ve[i], a = a, b = b, beta = beta, )
  #  oo[[i]] <- o
  #}

  pp <- unlist(lapply(q, function(x, ... ){
    se <-  seq(0, x, by = 1)
    px <- sum(dLKD(se, a=a, b = b, beta = beta), na.rm = T)
  }))

  if(log) pp <- log(pp)

  return(pp)
}



#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Probabilidades acumuladas da LKD.
#' @param x quantis de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param log se TRUE retorna o log das probabilidades
#' @param ... passagem de argumentos
#'
#' @details Para lidar com fatoriais envolvendo números decimais utilizamos a função Gamma, além disso, visando resolver problemas de números muito grandes nos fatoriais, usamos a aproximação de Stirling.
#'
#' @examples dLKD(x = 0:10, a = 5, b = 0.05, beta = 0.1)
#' @export
dLKD <- function(x, a = 5, b = 0.02, beta = 0.3, log = FALSE, ...){
  #if(a <= 0) stop("a deve ser > 0")
  #if(b <= -1) stop("b deve ser > -1")
  #if(p <= -b | p >= 1) stop("p deve ser > -b e menor que 1")
  if(a > 0 & b > -1 & (beta > -b & beta < 1)) {
    stirling <- function(z){
      sqrt(2*pi*z)*(z/exp(1))^z
    }

    o <- sapply(x, function(i){
      f1 <- a/beta + b*i/beta
      g1 <- gamma(f1 + i + 1)
      if(!is.finite(g1)){
        g1 <- stirling(f1 + i + 1)
        if(!is.finite(g1)) g1 <- .Machine$double.xmax
      }

      g2 <- gamma(i+1)
      if(!is.finite(g2)){
        g2 <- stirling(i+1)
        if(!is.finite(g2)) g2 <- .Machine$double.xmax
      }

      g3 <- gamma(f1 + 1)
      if(!is.finite(g3)){
        g3 <- stirling(f1 + 1)
        if(!is.finite(g3)) g3 <- .Machine$double.xmax
      }

      oo <- ((a/beta)/(f1 + i))*(g1/(g2*g3))*beta^(i)*(1-beta)^f1
      return(oo)
    })
    if(log) o <- log(o)
    return(o)
  } else {
    cat('indefinida nestes parametros!\n')
    return(NA)
  }
}


##-----------------------------------------------------------------------------
## Usando recursos númericos para simular da normal via método da
## transformação integral da probabilidade.

#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Probabilidades acumuladas da LKD.
#' @param p vetor de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param ... passagem de argumentos
#'
#' @details Gerar quantis da LKD utilizando a transformação integral de probabilidade através de aproximação numérica.
#'
#' @examples qLKD(10, a = 5, b = 0.05, beta = 0.1)
#' @export
qLKD <- function(p, a = 5, b = 0.02, beta = 0.3, ...){
  prec = 5*a
  xx <- 0:prec
  Fx <- pLKD(q = xx, a, b, beta)
  # aproxFx <- approxfun(x = xx, y = Fx)
  invFx <- approxfun(x = Fx, y = xx, rule = 2)
  #u <- dLKD(x = p, a = a, b = b, beta = beta)
  #cat("Fx", u, "\t")
  #u <- runif(n, min = 0, max = 1)
  rand <- as.integer(invFx(p))
  return(rand)
}


#' Distribuicao Katz Lagrangiana Generalizada (GLKD)
#'
#' Probabilidades acumuladas da GLKD.
#' @param x quantis de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param c parametro 3
#' @param beta parametro 4, tambem de nome beta
#' @param log se TRUE retorna o log das probabilidades
#' @param ... passagem de argumentos
#'
#' @details Para lidar com fatoriais envolvendo números decimais utilizamos a função Gamma, além disso, visando resolver problemas de números muito grandes nos fatoriais, usamos a aproximação de Stirling.
#'
#' @examples dGLKD(x = 0:10, a = 5, b = 0.05, c = 0.5, beta = 0.1)
## @export
dGLKD <- function(x, a = 5, b = 0.02, c = 0.5, beta = 0.3, log = FALSE, ...){
  #if(a <= 0) stop("a deve ser > 0")
  #if(b <= -1) stop("b deve ser > -1")
  #if(p <= -b | p >= 1) stop("p deve ser > -b e menor que 1")
  if(a > 0 & c > 0 & b > -c & (beta > 0 & beta < 1)) {
    stirling <- function(z){
      sqrt(2*pi*z)*(z/exp(1))^z
    }

    o <- sapply(x, function(i){
      f1 <- a/c + b*i/c
      g1 <- gamma(f1 + i + 1)
      if(!is.finite(g1)){
        g1 <- stirling(f1 + i + 1)
        if(!is.finite(g1)) g1 <- .Machine$double.xmax
      }

      g2 <- gamma(i+1)
      if(!is.finite(g2)){
        g2 <- stirling(i+1)
        if(!is.finite(g2)) g2 <- .Machine$double.xmax
      }

      g3 <- gamma(f1 + 1)
      if(!is.finite(g3)){
        g3 <- stirling(f1 + 1)
        if(!is.finite(g3)) g3 <- .Machine$double.xmax
      }

      oo <- ((a/c)/(f1 + i))*(g1/(g2*g3))*beta^(i)*(1-beta)^f1
      return(oo)
    })
    if(log) o <- log(o)
    return(o)
  } else {
    cat('indefinida nestes parametros!\n')
    return(NA)
  }
}

## Formas da LKD
## 1. Katz(a, beta)
## se b = 0

## 2. Binomial(n, theta)
## se 0 < beta = theta < 1; a = n*theta e b = -theta
## ou se b = 0; beta < 0; theta = beta*(beta-1)^(-1)


## 3. Poisson(a)
## se b = 0; beta --> 0 (beta tende a 0)

## 3. Poisson Generalizada(a, lambda), capítulo 9
## se b = 0; beta --> 0 (beta tende a 0)

## 5. Binomial Negativa BN(k, theta)
## se b = 0; 0 < beta = theta < 1; a/beta = k
## ou se b = theta(1-theta)^(-1); a = k*theta*(1-theta)^(-1)

## 6. Binomial Negativa Generalizada GNBD, capítulo 10
## se 0 < beta = theta < 1; a = n*theta e b = (m-1)*theta
## ou se beta < 0; beta = -theta(1-theta)^(-1); a = n*theta*(1-theta)^(-1) e b = m*theta*(1-theta)^(-1)

#' Dados os parametros verifica se é lkd
#'
#' @export
is.lkd <- function(n, a, b, beta, c = NULL){
	if((a > 0) & (b > -beta) & (beta < 1)) TRUE else FALSE
}

#' Dados os parametros verifica se a lkd é uma katz
#'
#' @export
is.kats <- function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (b == 0)) TRUE else FALSE
}

#' Dados os parametros verifica se a lkd é uma Binomial(N, p = beta)
#'
#' @export
is.bn <-  function(n, a, b, beta, c = NULL){
  theta <- beta
  if(beta > 0 & theta < 1 & a == n*theta & b == -theta) TRUE else FALSE
}

#' Dados os parametros verifica se a lkd é uma Poisson(theta = a)
#'
#' @export
is.po <-  function(n, a, b, beta, c = NULL){
	if((b == 0) & (round(beta, 10) == 0)) TRUE else FALSE
}

#' Dados os parametros verifica se a lkd é uma Binomial Negativa(N, p = beta)
#'
#' @export
is.nbn <-  function(n, a, b, beta, c = NULL){
	if((b == 0) & (beta > 0 & beta < 1)) TRUE else FALSE
}

#' Dados os parametros verifica se a lkd é uma Binomial Negativa Generalizada(N, beta, a)
#'
#' @export
is.gnbn <-  function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (beta > 0 & beta < 1) & (a == n*beta)) TRUE else FALSE
}


#' Estimativas de momento da lkd
#'
#' @export
est_mm <- function(y){
  xb <- mean(y)
  S2 <- var(y)
  S3 <- sum((y-xb)^3)/(length(y)-1)
  A <- (3*S2^2-S3*xb)^2 / (xb-S2^3)
  if(A < 4) {
    stop("A < 4, nao ha estimadores de momento para esta amostra.")
  } else {
    m_a <- 0.5*xb^(3/2)*(sqrt(A) + c(-1,1)*sqrt(A-4)) * S2^(-0.5)
    m_b <- -1 + 0.5*(sqrt(A) + c(-1,1)*sqrt(A-4))*(sqrt(A) - sqrt(xb/S2))
    m_beta <- 2-0.5*(A + c(-1,1)*sqrt(A*(A-4)))
  }
  #o <- c(m_a[1], m_b[1], m_beta[1])
  return(c(m_a = m_a[1], m_b = m_b[1], m_beta = m_beta[1]))
}

#' Estimativas com bese em momemento e zero-class para a lkd
#'
#' @export
est_mzcf <- function(y, a, beta, ...){
  f0 <- (1-beta)^(a/beta)
  #f0 <- sum(y==0)/ length(y)
  lf0 <- log(f0)

  # log de f0 não pode ser zero, aplicando o log do menor numero real posivel
  lf0[is.infinite(lf0)] <- log(.Machine$double.eps)

  xb <- mean(y, na.rm = T) #a/(1-b-beta)
  S2 <- var(y, na.rm = T) #a*(1-beta)/(1-b-beta)^3

  f2 <- function(x, ...){
    c((1-x)*(log(1-x))^2 -x^2*S2*(lf0)^2/xb^3)
  }

  beta_bar <- rootSolve::uniroot(f2, c(-.Machine$double.eps, 0.999999))$root

  a_bar <- (xb^3*(1-beta_bar)/S2)^(1/2)
  b_bar <- 1-beta_bar - a_bar/xb
  return(c(m_a = a_bar, m_b = b_bar, m_beta = beta_bar))
}

#' Estimativas de Maxima log-verossimilhança da lkd
#'
#' @export
est_mle <- function(y, a, beta, ...){
  n <- length(y)
  ni <- table(y)
  #f0 <- sum(y==0) / length(y)
  xb <- mean(y, na.rm = T) # a/(1-b-beta)
  #f0 <- (1-beta)^(a/beta)
  #f0 <- sum(y==0)/ length(y)
  #lf0 <- log(f0)

  # log de f0 não pode ser zero, aplicando o log do menor numero real posivel
  #lf0[is.infinite(lf0)] <- log(.Machine$double.eps)
  model <- function(par, ...){
    c(
      F1 = (n*xb*log(1-par[2]))/par[2] + sum(unlist(sapply(2:length(ni), function(i) {
        sapply(1:(i-1), function(j) {
          i*ni[i] / (xb*(1-par[1]-par[2]) + par[1]*i + par[2]*j)
        })
      })), na.rm = T)
      ,
      F2 = ((-n*xb*(1-par[2])*log(1-par[2]))/par[2]) - n*xb/par[2] + sum(unlist(sum(unlist(sapply(2:length(ni), function(i) {
        sapply(1:(i-1), function(j) {
          j*ni[i] / (xb*(1-par[1]-par[2]) + par[1]*i + par[2]*j)
        })
      })))), na.rm = T)
    )
  }

  if(sum(y==0) == 0){
    ## inicia com estimadores de momento caso não tenha zero-class
    o <- try(est_mm(y)[2:3])
    if(class(o)[1] != "try-error"){
      ini <- o
    } else {
      ## caso tenha zero-class, usa como start estimativas de momento e zero-class
      o <- try(est_mzcf(y)[2:3])
      if(class(o)[1] != "try-error"){
        ini <- o
      } else {
        ## caso conatrário usa um chute genérico para b e beta
        ini <- c(m_b = sd(y), m_beta = 0.5)
      }
    }
  }

  b_beta <- unname(rootSolve::multiroot(model, ini[2:3])$root)
  a <- xb*(1-b_beta[1] - b_beta[2])
  return(c(a = a, b = b_beta[1], beta = b_beta[2]))
}


#' Estimativas de Maxima log-verossimilhança da lkd em outra abordagem
#'
#' @export
est_mle_v2 <- function(y){
  n <- length(y)
  ni <- table(y)
  #f0 <- sum(y==0) / length(y)
  xb <- mean(y, na.rm = T) # a/(1-b-beta)
  start <- est_mzcf(y)

  model <- function(pars, ...){
    b <- pars[1]; beta <- pars[2]
    c(
      F1 =  c(n * xb * (1 - beta) / beta + sum(unlist(sapply(2:length(ni), function(i) {
        i * ni[i] * digamma(i - (xb * b + xb * beta - b * i - xb) / b) / b - i * ni[i] * digamma(1 - (xb * b + xb * beta - b * i - xb) / b) / b
      })), na.rm = TRUE))

      ,
      F2 = c(-n * xb * (1 - beta) * log1p(-beta) - n * xb / beta + sum(unlist(sapply(2:length(ni), function(i) {
        i * ni[i] * digamma(i - (xb * b + xb * beta - b * i - xb) / b) / b - i * ni[i] * digamma(1 - (xb * b + xb * beta - b * i - xb) / b) / b
      })), na.rm = TRUE))
    )
  }

  b_beta <- unname(rootSolve::multiroot(model, c(start[2], start[3]))$root)
  a <- xb*(1-b_beta[1] - b_beta[2])
  return(c(a = a, b = b_beta[1], beta = b_beta[2]))
}

#' Função gradiente analítica da lkd(a,b,beta)
#'
#' @export
grad_lkd <- function(pars, ...){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)

  ## Computando dy/a
  a1 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * digamma(i + (b * i + a) / beta) / beta - ni[i] * digamma(1 + (b * i + a) / beta) / beta
  })))
  xa <- -(n - n0) / a ^ 2 * log(a) + (n - n0) / a ^ 2 + n / beta * log1p(-beta) + a1

  ## Computando dy/b
  b1 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * digamma(i + (b * i + a) / beta) * i / beta - ni[i] * digamma(1 + (b * i + a) / beta) * i / beta
  })))
  xb <- n * xbar / beta * log1p(-beta) + b1

  ## Computando dy/beta

  s1 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * (digamma((b * i + i * beta + a) / beta) * b * i - digamma((b * i + a + beta) / beta) * b * i + digamma((b * i + i * beta + a) / beta) * a - digamma((b * i + a + beta) / beta) * a - i * beta + beta)
  })))

  s2 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * (digamma((b * i + i * beta + a) / beta) * b * i - digamma((b * i + a + beta) / beta) * b * i + digamma((b * i + i * beta + a) / beta) * a - digamma((b * i + a + beta) / beta) * a - i * beta + beta)
  })))

  xbeta <- -(log1p(-beta) * xbar * b * beta * n - log1p(-beta) * xbar * b * n + log1p(-beta) * a * beta * n - xbar * b * beta * n - log1p(-beta) * a * n - a * beta * n + s1 * beta - s2) / beta ^ 2 / (-1 + beta)

  return(c(xa, xb, xbeta))
}
grad_lkd <- compiler::cmpfun(grad_lkd)

#' Função Hessiana analítica da lkd(a,b,beta)
#'
#' @export
hessian_lkd <- function(pars, ...){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)

  ## Obtendo dy^2/a^2
  s1 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })))

  aa <- 2 * (n - n0) / a ^ 3 * log(a) - 3 * (n - n0) / a ^ 3 + s1

  ## Obtendo dy^2/a*b
  ab <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })))

  ## Obtendo dy^2/a*beta
  s2 <- sum(unlist(sapply(2:length(ni), function(i) {
    -ni[i] * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) - ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) + ni[i] * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })))

  abeta <- -n / beta ^ 2 * log1p(-beta) - n / beta / (1 - beta) + s2

  ## Obtendo dy^2/b*a
  s3 <- sum(unlist(sapply(2:length(ni), function(i) {
    -ni[i] * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) - ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) + ni[i] * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })))

  abeta <- -n / beta ^ 2 * log1p(-beta) - n / beta / (1 - beta) + s3

  ## Obtendo dy^2/b*a
  ba <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * i / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] * i / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })))

  ## Obtendo dy^2/b*b
  bb <-  sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * i ^ 2 / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] * i ^ 2 / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })))

  ## Obtendo dy^2/b*beta
  s4 <- sum(unlist(sapply(2:length(ni), function(i) {
    -ni[i] * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) * i - ni[i] * digamma(i + (b * i + a) / beta) * i / beta ^ 2 + ni[i] * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) * i + ni[i] * digamma(1 + (b * i + a) / beta) * i / beta ^ 2
  })))

  bbeta <- -n * xbar / beta ^ 2 * log1p(-beta) - n * xbar / beta / (1 - beta) + s4

  ## Obtendo dy^2/betaa
  s5 <- sum(unlist(sapply(2:length(ni), function(i) {
    -ni[i] * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) - ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) + ni[i] * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })))

  betaa <- -n / beta ^ 2 * log1p(-beta) - n / beta / (1 - beta) + s5

  ## Obtendo dy^2/betab
  s6 <- sum(unlist(sapply(2:length(ni), function(i) {
    -ni[i] * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) * i - ni[i] * digamma(i + (b * i + a) / beta) * i / beta ^ 2 + ni[i] * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) * i + ni[i] * digamma(1 + (b * i + a) / beta) * i / beta ^ 2
  })))

  betab <- -n * xbar / beta ^ 2 * log1p(-beta) - n * xbar / beta / (1 - beta) + s6

  ## Obtendo dy^2/betabeta
  s7 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * (-i / beta ^ 2 + (b * i + a) ^ 2 / beta ^ 4 * trigamma(i + (b * i + a) / beta) + 2 * digamma(i + (b * i + a) / beta) * (b * i + a) / beta ^ 3) - ni[i] * (-1 / beta ^ 2 + (b * i + a) ^ 2 / beta ^ 4 * trigamma(1 + (b * i + a) / beta) + 2 * digamma(1 + (b * i + a) / beta) * (b * i + a) / beta ^ 3)
  })))

  betabeta <- 2 * n * (xbar * b + a) / beta ^ 3 * log1p(-beta) + 2 * n * (xbar * b + a) / beta ^ 2 / (1 - beta) - n * (xbar * b + a) / beta / (1 - beta) ^ 2 + s7

  o <- matrix(c(aa, ab, abeta,
    ba, bb, bbeta,
    betaa, betab, betabeta),nrow=3,ncol=3, byrow = TRUE)
  return(o)
}
hessian_lkd <- compiler::cmpfun(hessian_lkd)

#' Função de log-verossimilhança altenativa da lkd(a,b,beta)
#'
#' @export
likeli_lkd <- function(pars, y){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)

  s0 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * log(factorial(i))
  })), na.rm = TRUE)

  s1 <- sum(unlist(sapply(2:length(ni), function(i) {
    ni[i] * (i * log(beta) + log(gamma(i + (b * i + a) / beta))) - ni[i] * (log(beta) + log(gamma(1 + (b * i + a) / beta)))
  })), na.rm = TRUE)

  ll <- (n - n0) / a * log(a) + n * (xbar * b + a) / beta * log1p(-beta) - s0 + s1
  return(ll)
}
likeli_lkd <- compiler::cmpfun(likeli_lkd)
