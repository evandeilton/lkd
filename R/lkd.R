## Exports
# @importFrom rootSolve uniroot uniroot.all multiroot
# @importFrom stats var sd
# @export
#stats::sd
# @export
#stats::sum
# @export
#stats::var
#'@importFrom stats approxfun dbinom dnbinom dpois pbinom pnbinom ppois qbinom qnbinom qpois rbinom rnbinom rpois runif sd uniroot var
# @importFrom MASS ginv
# @import rootSolve
NULL

#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Esta função retorna probabilidades da LKD. Podendo ser: lkd, kd, binomial, binomial negativa, binomial negativa generalizada, Poisson e Poisson Generalizada de acordo com os parãmetros a, b e beta.
#'
#' @param x um vetor de observações de uma lkd(a, b, beta)
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param size utlizado quando a lkd deriva Binomial ou Binomial Negativa
#' @param f_bino se TRUE retorna dlkd na forma fatorial, se FALSE retorna a forma função gamma.
#' @param log se TRUE retorna o log das probabilidades
#' @param ... passagem de argumentos
#'
#' @details Para lidar com fatoriais envolvendo números decimais utilizamos a função Gamma. Além disso, todas as contas são feitas em log e ao final retornam em exponencial.
#'
#' @return vetor de probabilidades da lkd(a, b, beta)
#'
#' @examples
#' X <- 0:20
#' sum(dlkd(X, 2, 0.1, 0.2, log = FALSE))
#' sum(dlkd(X, 2, 0, 0.2))
#' sum(dlkd(X, 2, 0.7, 0.9, log = FALSE, f_bino = FALSE))
#' sum(dlkd(X, 2, 0, 0, log = FALSE))
#' plot(function(x) dlkd(x, 20, 0.5, 0.9), 0, 50, main = "LKD", type = "h")
#'
#' @export
dlkd <- function(x, a = 5, b = 0.02, beta = 0.3, size = 1, f_bino = FALSE, log = FALSE, ...){
  if(a < 0) stop("a deve ser > 0")
  if(beta > 1) stop("beta deve ser < 1")
  if(b < -beta) stop("b > que -beta")

  X <- as.numeric(x)
  n <- length(X)

  if(is.bn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial(n, size, beta).
      Retornando probabilidades Binomial(n, size = 1, prob = beta).
      Edite o size para outras probabilidades Binomiais! \n")
    XX <- dbinom(x, size, beta, log = log)
  } else if (is.po(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Poisson(lambda = a)
      Retornando amostras Poisson \n")
    XX <- dpois(x, a, log = log)
  } else if (is.nbn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial Negativa BN(n, size, beta).
      Retornando probabilidades BN(n, size = 1, prob = beta).
      Edite o size para mudar o numero de sucessos e obter outras probabilidades! \n")
    XX <- dnbinom(x, size, beta, log = log)
  } else {
    ## Se a forma da lkd for binomial, senão usa gama para computar fatoriais
    if(f_bino){
      XX <- (log(choose(((b + beta) * x + a) / beta, x)) * beta - log((b + beta) * x + a) * beta + (b * x + a) * log1p(-beta) + beta * (x * log(beta) + log(a))) / beta
      if(!log) {
        XX <- exp(XX)
      }
    } else {
      XX <- (lgamma(((b + beta) * x + a) / beta) * beta - lgamma((b * x + a + beta) / beta) * beta + (b * x + a) * log1p(-beta) + (-lgamma(x + 1) + (x - 1) * log(beta) + log(a)) * beta) / beta
      if(!log) {
        XX <- exp(XX)
      }
    }
  }
  return(XX)
}


#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Retorna probabilidades acumuladas de uma lkd(a, b, beta)
#'
#' @param q vetor de quantis de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param size utlizado quando a lkd deriva Binomial ou Binomial Negativa
#' @param log se TRUE retorna o log das probabilidades
#' @param f_bino se TRUE retorna dlkd na forma fatorial, se FALSE retorna a forma função gamma.
#' @param ... passagem de argumentos
#' @details todo
#'
#' @return vetor de probabilidades acumuladas para cada quartil da lkd(a, b, beta)
#'
#' @examples
#' plkd(0:5, 2, -0.2, 0.2)
#' plkd(0:5, 2, -0.2, 0.2, log = TRUE)
#' plkd(0:5, 2, -0.2, 0.2, log = TRUE, f_bino = FALSE)
#' plkd(0:5, 2, 0, 0.2)
#' plkd(0:5, 2, 0.7, .9)
#' plot(function(x) plkd(x, 20, 0.5, 0.9), 0, 50, main = "LKD", type = "h")
#'
#' @export
plkd <- function(q, a = 5, b = 0.02, beta = 0.3, size = 1, f_bino = FALSE, log = FALSE, ...){
  if(a < 0) stop("a deve ser > 0")
  if(beta > 1) stop("beta deve ser < 1")
  if(b < -beta) stop("b > que -beta")
  n <- length(q)
  ## Se for binomial,toma amostra da binomial(n, tehe)
  if(is.bn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial(n, size, beta).
      Retornando probabilidades Binomial(n, size = 1, prob = beta).
      Edite o size para outras probabilidades Binomiais! \n")
    Q <- pbinom(q, size, beta, log.p = log)
  } else if (is.po(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Poisson(lambda = a)
      Retornando amostras Poisson \n")
    Q <- ppois(q, a, log.p = log)
  } else if (is.nbn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial Negativa BN(n, size, beta).
      Retornando probabilidades BN(n, size = 1, prob = beta).
      Edite o size para mudar o numero de sucessos e obter outras probabilidades! \n")
    Q <- pnbinom(q, size, beta, log.p = log)
  } else {
    Q <- unlist(lapply(q, function(x, ... ){
      se <-  seq(0, x, by = 1)
      px <- sum(dlkd(se, a=a, b = b, beta = beta, log = log, f_bino = f_bino), na.rm = TRUE)
    }))
  }
  return(Q)
}


#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Geração de números aleatórios da LKD pelo método da inversa da função de probabilidade acumulada.
#' Ocasionalmente a lkd pode derivar outras ditribuições. Por hora temos: lkd, kd, binomial, binomial negativa e Poisson. Outras serão configuradas em outras versões. No caso de recair em uma distribuição conhecida, o algoritimo de geração de numeros aleatórios destas distribuições é utilizado. Vide \code{\link{rbinom}}, \code{\link{rnbinom}} e \code{\link{rpois}}. Demais casos recam no algoritimo desenvolvido com base em Consul e Famoye (2006).
#'
#' @param n número de observacções desejadas da lkd(a, b, beta)
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3
#' @param size utlizado quando a lkd deriva Binomial ou Binomial Negativa
#' @param log se TRUE retorna o log das probabilidades
#' @param acum se TRUE retorna a probabilidade acumulada no quatil X, senão retorna X como uma realização da lkd(a, b, beta)
#' @param ... passagem de argumentos
#'
#' @details Esta função utiliza uma função interna que determina amostras da lkd(a, b, beta) pelo método da inversa da função de probabilidade acumulada.
#'
#' @return vetor de observações de uma lkd(a, b, beta)
#'
#' @examples
#' # rlkd(10, 2, -0.2, 0.2)
#' # rlkd(10, 2, -0.2, 0.3)
#' # rlkd(100, 2, 0, 0.2)
#' # rlkd(100, 2, 0.7, 6)
#' # x <- rlkd(500, 2, 0.2, 0.3)
#' # hist(x, main = 'Amostras da lkd(2, 0.2, 0.3)')
#' @export
#'
rlkd <- function(n, a = 5, b = 0.02, beta = 0.3, size = 1, log = FALSE, acum = FALSE, ...){
  if(a < 0) stop("a deve ser > 0")
  if(beta > 1) stop("beta deve ser < 1")
  if(b < -beta) stop("b > que -beta")

  ## Se for binomial,toma amostra da binomial(n, beta)
  if(is.bn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial(n, size, beta).
      Retornandoamostras Binomial(n, size = 1, prob = beta).
      Edite o size para amostras Binomiais! \n")
    S <- rbinom(n, size, beta)
  } else if (is.po(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Poisson(lambda = a)
      Retornando amostras Poisson \n")
    S <- rpois(n, a)
  } else if (is.nbn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial Negativa BN(n, size, beta).
      Retornando amostras BN(n, size = 1, prob = beta).
      Edite o size para mudar o numero de sucessos! \n")
    S <- rnbinom(n, size, beta)
  } else {
    ve <- sort(runif(n, min = 0, max = 1))
    oo <- c()
    for(i in seq_along(ve)){
      o <- .xlkd(u = ve[i], a = a, b = b, beta = beta)
      oo[[i]] <- o
    }
    S <- sample(oo, replace = FALSE)
  }
  if(log) {
    S <- log(S)
  }
  return(S)
}

## @export
.xlkd <- function(u, a, b, beta, acum = FALSE){
  if(a < 0) stop("a deve ser > 0")
  if(beta > 1) stop("beta deve ser < 1")
  if(b < -beta) stop("b > que -beta")

  X <- 0
  P0 <- log((1-beta)^(a/beta))
  P0 <- exp(P0)
  S <- P <- P0
  U <- u
  if(U <= S){
    if(acum) return(S) else return(X)
  }
  X <- X + 1
  P <- log(a*(1-beta)^(b/beta)*S)
  P <- exp(P)
  S <- S + P
  if(U <= S){
    if(acum) return(S) else return(X)
  }

  while(U > S){
    X <- X + 1
    c1 <- log((a + b*(X+1) + beta*X) / (X + 1))
    c1 <- exp(c1)

    c2 <- log((1-beta)^(b/beta))
    c2 <- exp(c2)

    c3 <- log(1 + (b/(a + b*X + beta*(1:(X-1)))))
    c3 <- sum(exp(c3), na.rm = TRUE)

    P <- c1*c2*c3*S
    S <- S + P
    #cat("x: ", X, "u: ", U, "P: ", P, "s: ", S, "\n")
  }
  if(acum) return(S) else return(X)
}


#' Distribuicao Katz Lagrangiana (LKD)
#'
#' Inversa das probabilidades acumuladas da lkd(a, b, beta) obtida por simulação via método da transformação integral da probabilidade.
#'
#' @param p vetor de probabilidade da LKD
#' @param a parametro 1
#' @param b parametro 2
#' @param beta parametro 3, tambem de nome beta
#' @param ... passagem de argumentos
#'
#' @details Gerar quantis da LKD utilizando a transformação integral de probabilidade através de aproximação numérica.
#'
#' @examples
#' p <- plkd(0:10, 2, 0.1, 0.2, log = FALSE)
#' qlkd(p, 2, 0.1, 0.2, log = FALSE)
#' qlkd(ppois(0:5, lambda = 2), 2, 0, 0, log = FALSE)
#' par(mfrow = c(1,2))
#' plot(function(x) plkd(x, 2, 0.1, 0.2), ylim = c(0,1), xlim = c(0, 15),
#'  main = "Acumulada LKD", type = "l")
#' plot(function(x) qlkd(x, 2, 0.1, 0.2), ylim = c(0,15), xlim = c(0, 1),
#'  main = "Inversa Acumulada LKD", type = "l")
#' par(mfrow = c(1,1))
#'
#' @export
qlkd <- function(p, a = 5, b = 0.02, beta = 0.3, log = FALSE, ...){
  if(a < 0) stop("a deve ser > 0")
  if(beta > 1) stop("beta deve ser < 1")
  if(b < -beta) stop("b > que -beta")
  n <- length(p)

  ## Se for binomial,toma amostra da binomial(n, beta)
  if(is.bn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial(n, size, beta).
      Retornandoamostras Binomial(n, size = 1, prob = beta).
      Edite o size para amostras Binomiais! \n")
    Q <- qbinom(p, size, beta, log.p = log)
  } else if (is.po(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Poisson(lambda = a)
      Retornando amostras Poisson \n")
    Q <- qpois(p, a, log.p = log)
  } else if (is.nbn(n, a, b, beta)){
    cat("Log: Esta lkd derivou uma Binomial Negativa BN(n, size, beta).
      Retornando amostras BN(n, size = 1, prob = beta).
      Edite o size para mudar o numero de sucessos! \n")
    Q <- qnbinom(p, size, beta, log.p = log)
  } else {
    prec = 50*a
    xx <- 0:prec
    Fx <- plkd(q = xx, a, b, beta)
    # aproxFx <- approxfun(x = xx, y = Fx)
    invFx <- approxfun(x = Fx, y = xx, rule = 2)
    #u <- dlkd(x = p, a = a, b = b, beta = beta)
    #cat("Fx", u, "\t")
    #u <- runif(n, min = 0, max = 1)
    Q <- as.integer(invFx(p))
    if(log) Q <- log(Q)
  }
  return(Q)
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
lkd_est_mm <- function(y){
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
#' @export
lkd_est_mzcf <- function(y, a, beta, ...){
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

  beta_bar <- uniroot(f2, c(-.Machine$double.eps, 0.999999))$root

  a_bar <- (xb^3*(1-beta_bar)/S2)^(1/2)
  b_bar <- 1-beta_bar - a_bar/xb
  return(c(m_a = a_bar, m_b = b_bar, m_beta = beta_bar))
}

#' Estimativas de Maxima log-verossimilhança da lkd
#'
#' @export
lkd_est_mle <- function(y, a, beta, ...){
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
    o <- try(lkd_est_mm(y)[2:3])
    if(class(o)[1] != "try-error"){
      ini <- o
    } else {
      ## caso tenha zero-class, usa como start estimativas de momento e zero-class
      o <- try(lkd_est_mzcf(y, a, beta)[2:3])
      if(class(o)[1] != "try-error"){
        ini <- o
      } else {
        ## caso conatrário usa um chute genérico para b e beta
        ini <- c(m_b = sd(y), m_beta = 0.5)
      }
    }
  } else {
    ini <- try(lkd_est_mzcf(y, a, beta)[2:3])
  }

  b_beta <- unname(rootSolve::multiroot(model, ini)$root)
  a <- xb*(1-b_beta[1] - b_beta[2])
  return(c(a = a, b = b_beta[1], beta = b_beta[2]))
}


#################################### abordagens com CAS ##################################

################### Funções para log-verossimilhança de Consul t Famoye ###################

lkd_log_vero_consul <- function(pars, ...){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  y <- y
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)

  s1 <- sapply(2:length(ni), function(i) {
    ni[i] * log(factorial(i))
  })

  #ii <- !is.finite(s1)
  #s1[ii] <- log(.Machine$double.xmin * (1.15e-16))
  s1 <- sum(unlist(s1), na.rm = T)

  s2 <- sapply(2:length(ni), function(i) {
    sapply(1:(i-1), function(j) {
      ni[i] * log(b * i + beta * j + a)
    })
  })
  s2 <- sum(unlist(s2), na.rm = T)

  o <- (n - n0) * log(a) + n * (b * xbar + a) / beta * log1p(-beta) - s1 + s2

  return(o)
}
#' Função de log-verosimilhança
#'
#' @description calcula a função de log-verossimilhança presente em Consul e Famoye (2016).
#' Esta função é pré-compilada para tornar sua execução mais rápida. Ela consome um vetor
#' de observações de uma lkd(a, b, beta) de nome y que deve estar no ambiente do R.
#'
#' @param pars vetor de parametros
#' @param ... passagem de argumentos
#' @export
lkd_log_vero_consul <- compiler::cmpfun(lkd_log_vero_consul)


lkd_gradient_log_vero_consul <- function(pars, ...){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)
  y <- y

  ## obtendo dy/a
  s1 <- sapply(2:length(ni), function(i) {
    ni[i] * digamma(i + (b * i + a) / beta) / beta - ni[i] * digamma(1 + (b * i + a) / beta) / beta
  })
  s1 <- sum(unlist(s1), na.rm = T)
  aa <- (n - n0) / a + n / beta * log1p(-beta) + s1

  ## obtendo dy/b
  s2 <- sapply(2:length(ni), function(i) {
    ni[i] * i * digamma(i + (b * i + a) / beta) / beta - ni[i] * i * digamma(1 + (b * i + a) / beta) / beta
  })
  s2 <- sum(unlist(s2), na.rm = T)
  bb <- n * xbar / beta * log1p(-beta) + s2

  ## obtendo dy/beta
  s4 <- sapply(2:length(ni), function(i) {
    ni[i] / beta * i - (b * i + a) * ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 - 1 / beta * ni[i] + (b * i + a) * ni[i] * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })
  s4 <- sum(unlist(s4), na.rm = T)

  betabeta <- -n * (b * xbar + a) / beta ^ 2 * log1p(-beta) - n * (b * xbar + a) / beta / (1 - beta) + s4

  o <- c(aa, bb, betabeta)
  return(o)
}

#' Função Gradiente analítica
#'
#' @description calcula analiticamente o gradiente da função de log-verossimilhança presente em Consul e Famoye (2016). Por definição, a função gradiente é formada por todas a as derivadas analílicas de primeira ordem da função de log-verossimilhança.  Esta função é pré-compilada para tornar sua execução mais rápida. Ela consome um vetor de observações de uma lkd(a, b, beta) de nome y que deve estar no ambiente do R.
#' @return vetor das tres derivadas de primeira ordem
#' @param pars vetor de parametros
#' @param ... passagem de argumentos
#' @export
lkd_gradient_log_vero_consul <- compiler::cmpfun(lkd_gradient_log_vero_consul)

lkd_hessian_log_vero_consul <- function(pars, ...){
  a <- pars[1];b<-pars[2]; beta <- pars[3]
  xbar <- mean(y)
  n0 <- sum(y == 0)
  n  <- length(y)
  ni <- table(y)
  y <- y

  ## Obtendo d²y/a²
  s0 <- sapply(2:length(ni), function(i) {
    ni[i] / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })
  s0 <- sum(unlist(s0), na.rm = T)
  aa <- -(n - n0) / a ^ 2 + s0

  ## Obtendo d²y/ab
  s1 <- sapply(2:length(ni), function(i) {
    ni[i] / beta ^ 2 * i * trigamma(i + (b * i + a) / beta) - ni[i] / beta ^ 2 * i * trigamma(1 + (b * i + a) / beta)
  })
  ab <- sum(unlist(s1), na.rm = T)

  ## Obtendo d²y/abeta
  s2 <- sapply(2:length(ni), function(i) {
    -ni[i] / beta ^ 3 * (b * i + a) * trigamma(i + (b * i + a) / beta) - ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] / beta ^ 3 * (b * i + a) * trigamma(1 + (b * i + a) / beta) + ni[i] / beta ^ 2 * digamma(1 + (b * i + a) / beta)
  })
  s2 <- sum(unlist(s2), na.rm = T)
  abeta <- -n / beta ^ 2 * log1p(-beta) - n / beta / (1 - beta) + s2

  ## Obtendo d²y/ba
  s3 <- sapply(2:length(ni), function(i) {
    ni[i] / beta ^ 2 * i * trigamma(i + (b * i + a) / beta) - ni[i] / beta ^ 2 * i * trigamma(1 + (b * i + a) / beta)
  })
  ba <- sum(unlist(s3), na.rm = T)

  ## Obtendo d²y/b²
  s4 <- sapply(2:length(ni), function(i) {
    ni[i] * i ^ 2 / beta ^ 2 * trigamma(i + (b * i + a) / beta) - ni[i] * i ^ 2 / beta ^ 2 * trigamma(1 + (b * i + a) / beta)
  })
  bb <- sum(unlist(s4), na.rm = T)

  ## Obtendo d²y/bbeta
  s5 <- sapply(2:length(ni), function(i) {
    -ni[i] * i * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) - ni[i] * i * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] * i * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) + ni[i] * i * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })
  s5 <- sum(unlist(s5), na.rm = T)
  bbeta <- -n * xbar / beta ^ 2 * log1p(-beta) - n * xbar / beta / (1 - beta) + s5

  ## Obtendo d²y/betaa
  s7 <- sapply(2:length(ni), function(i) {
    -ni[i] / beta ^ 3 * (b * i + a) * trigamma(i + (b * i + a) / beta) - ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] / beta ^ 3 * (b * i + a) * trigamma(1 + (b * i + a) / beta) + ni[i] / beta ^ 2 * digamma(1 + (b * i + a) / beta)
  })
  s7 <- sum(unlist(s7), na.rm = T)
  betaa <- -n / beta ^ 2 * log1p(-beta) - n / beta / (1 - beta) + s7

  ## Obtendo d²y/betab
  s8 <- sapply(2:length(ni), function(i) {
    -ni[i] * i * (b * i + a) / beta ^ 3 * trigamma(i + (b * i + a) / beta) - ni[i] * i * digamma(i + (b * i + a) / beta) / beta ^ 2 + ni[i] * i * (b * i + a) / beta ^ 3 * trigamma(1 + (b * i + a) / beta) + ni[i] * i * digamma(1 + (b * i + a) / beta) / beta ^ 2
  })
  s8 <- sum(unlist(s8), na.rm = T)
  betab <-   -n * xbar / beta ^ 2 * log1p(-beta) - n * xbar / beta / (1 - beta) + s8

  ## Obtendo d²y/beta²
  s10 <- sapply(2:length(ni), function(i) {
    -ni[i] / beta ^ 2 * i + (b * i + a) ^ 2 * ni[i] / beta ^ 4 * trigamma(i + (b * i + a) / beta) + 2 * (b * i + a) * ni[i] * digamma(i + (b * i + a) / beta) / beta ^ 3 + 1 / beta ^ 2 * ni[i] - (b * i + a) ^ 2 * ni[i] / beta ^ 4 * trigamma(1 + (b * i + a) / beta) - 2 * (b * i + a) * ni[i] * digamma(1 + (b * i + a) / beta) / beta ^ 3
  })
  s10 <- sum(unlist(s10), na.rm = T)
  betabeta <-   2 * n * (b * xbar + a) / beta ^ 3 * log1p(-beta) + 2 * n * (b * xbar + a) / beta ^ 2 / (1 - beta) - n * (b * xbar + a) / beta / (1 - beta) ^ 2 + s10

  o <- matrix(c(aa, ab, abeta, ba, bb, bbeta, betaa, betab, betabeta), nrow=3,ncol=3, byrow = TRUE)

  return(o)
}

#' Função Hessiana analítica
#'
#' @description calcula analiticamente o Hessiano da função de log-verossimilhança presente em Consul e Famoye (2016). Por definição, a matriz Hessiana é a matrix das derivadas analílicas de segunda ordem da função de log-verossimilhança.  Esta função é pré-compilada para tornar sua execução mais rápida. Ela consome um vetor de observações de uma lkd(a, b, beta) de nome y que deve estar no ambiente do R.
#' @param pars vetor de parametros
#' @param ... passagem de argumentos
#' @return matriz 3 x 3.
#' @export
lkd_hessian_log_vero_consul <- compiler::cmpfun(lkd_hessian_log_vero_consul)


#' Algoritmo Newton-Rapson matricial
#'
#' @description Função generica para obter estimativas por Newton-Rapson usando o gradiente e o hessiano da função. Se o sistema dormado pelo Hessiano for singular, ou seja, sem solução. Usa-se a inversa de Moore-Pen Rose como aproximação para as estimativas.
#'
#' @param initial valores do chute inicial
#' @param escore função escore da log-verossimilhança
#' @param hessiano função Hessiana da log-verossimilhança
#' @param tol tolerãncia para o erro entre uma estimativa e a anterior
#' @param max.iter numero de iterações alogortmo antes da convergência
#' @param n.dim número de dimensões (parãmetros) da log-verossimilhança
#' @param trace se TRUE esibe as estimativas em cada iteração
#' @param ... outros parametros
#' @return lista contendo dois elementos: o Hessiano da ultima iteração e as estimativas
#' @examples

#' # N <- 500
#' # y <- 0:(N-1)
#' # a <- 10; b <- 0.5; beta <- 0.05
#' # pars <- c(a, b, beta)
#' # y <- rlkd(N, a, b, beta)
#' # hist(y)
#' # start <- unname(lkd_est_mzcf(y, a, beta))
#' # k <- newton_raphson(initial = start
#' # , escore = lkd_gradient_log_vero_consul
#' # , hessiano = lkd_hessian_log_vero_consul
#' # , max.iter = 5000, n.dim = 3, tol = 0.00001)
#' @export
newton_raphson <- function(initial, escore, hessiano, tol=0.0001, max.iter = 500, n.dim = 3, trace = TRUE, ...){
  solucao <- matrix(NA, max.iter,n.dim)
  solucao[1,] <- initial
  for(i in 2:max.iter){
    HSS <- hessiano(initial)
    ## ESC mudou para as.numeric
    ESC <- as.numeric(t(escore(initial)))

    ## Adicionado invesa generalizada quando o sistema for singular
    ss <- try(solve(HSS,ESC))
    if(class(ss) != "try-error") {
      S <- ss
      solucao[i,] <- initial - S
    } else {
      S <- MASS::ginv(HSS) %*% ESC
      solucao[i,] <- initial - S
    }

    initial <- solucao[i,]
    tolera <- abs(solucao[i,] - solucao[i-1,])
    if(trace) cat("par =", solucao[i,], "\n")
    if(all(tolera < tol) == TRUE) break
  }
  saida <- list()
  saida[1][[1]] <- HSS
  saida[2][[1]] <- initial
  saida <<- initial
  return(saida)
}


