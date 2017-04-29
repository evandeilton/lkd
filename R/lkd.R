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

is.lkd <- function(n, a, b, beta, c = NULL){
	if((a > 0) & (b > -beta) & (beta < 1)) TRUE else FALSE
}

is.kats <- function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (b == 0)) TRUE else FALSE
}

is.bn <-  function(n, a, b, beta, c = NULL){
  theta <- beta
  if(beta > 0 & theta < 1 & a == n*theta & b == -theta) TRUE else FALSE
}

is.po <-  function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (b == 0) & (round(beta, 5) == 0)) TRUE else FALSE
}

is.nbn <-  function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (b == 0) & (beta > 0 & beta < 1)) TRUE else FALSE
}

is.gnbn <-  function(n, a, b, beta, c = NULL){
	if(is.lkd(n, a, b, beta, c) & (beta > 0 & beta < 1) & (a == n*beta)) TRUE else FALSE
}

