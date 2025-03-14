
### By J.A. Fiorucci
## juntanto os scripts obtidos do github do aluno: https://github.com/JG2108/ULARMA

#require("expint")


################################################################################
#### Erro Padrao Heissiano
SEfromHessian <- function(a, hessian=FALSE, silent=FALSE) {
  namesp <- colnames(a)
  mathessian <- a

  mathessian <- ifelse(mathessian==-Inf, -1E9, mathessian)
  mathessian <- ifelse(mathessian==+Inf, 1E9, mathessian)

  sigma  <- try(solve(mathessian), silent=TRUE)
  # Add all: 22/4/2020; any !
  if (inherits(sigma, "try-error")) {
    if (!silent) warning("Error in Hessian matrix inversion")
    mathessianx <- try(as.matrix(getFromNamespace("nearPD", ns="Matrix")(mathessian)$mat), silent=TRUE)
    # 29/1/2021
    if (inherits(mathessianx, "try-error")) {
      if (!silent) warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
      sigma  <- try(ginv(mathessian), silent=TRUE)
    } else {
      if (!silent) warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
      sigma  <- try(solve(mathessianx), silent=TRUE)
    }
  }

  # Add all: 22/4/2020
  if (!inherits(sigma, "try-error")) {

    if (all(diag(sigma)>=0)) {
      # méthode classique
      res <- sqrt(diag(sigma))

    } else {
      # Autre essai
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      res <- structure((matrix(rep(1, nrow(R)), nrow = 1, byrow = TRUE) %*% R)[1, ], .Names=colnames(mathessian))

      # Si j'ai des SE négatif... pas bon signe
      if (any(res<0)) {
        d <- diag(as.matrix(getFromNamespace("nearPD", ns="Matrix")(sigma)$mat))
        names(d) <- colnames(mathessian)
        res <- ifelse(d<0, NA, sqrt(d))
      }
      if (any(is.na(res))) {


        a <- sigma

        n = dim(a)[1];
        root = matrix(0,n,n);

        for (i in 1:n){
          sum = 0;
          if (i>1){
            sum = sum(root[i,1:(i-1)]^2);
          }

          x = a[i,i] - sum;

          if (x<0){
            x = 0;
          }

          root[i,i] = sqrt(x);

          if (i < n){
            for (j in (i+1):n){

              if (root[i,i] == 0){
                x=0;
              }
              else{
                sum = 0;
                if (i>1) {
                  sum = root[i,1:(i-1)] %*% t(t(root[j,1:(i-1)]))
                }
                x = (a[i,j] - sum)/root[i,i];
              }

              root[j,i] = x;
            }
          }
        }
        colnames(root) <- rownames(root) <- colnames(mathessian)

        pseudoV <- root %*% t(root)

        d <- diag(pseudoV)

        if (any(d != 0) & all(d >= 0)) {
          res <- sqrt(d)

          if (!silent) warning("Estimates using pseudo-variance based on Gill & King (2004)")

        } else {
          if (!silent) warning("Approximation of Cholesky matrix based on Rebonato and Jackel (2000)")
          if (!silent) warning("Estimates using pseudo-variance based on Gill & King (2004)")

          # The paper by Rebonato and Jackel, “The most general methodology for creating a valid correlation matrix for risk management and option pricing purposes”, Journal of Risk, Vol 2, No 2, 2000, presents a methodology to create a positive definite matrix out of a non-positive definite matrix.
          # FP Brissette, M Khalili, R Leconte, Journal of Hydrology, 2007, “Efficient stochastic generation of multi-site synthetic precipitation data”
          # GA Baigorria, JW Jones, Journal of Climate, 2010, “GiST: A stochastic model for generating spatially and temporally correlated daily rainfall data”
          # M Mhanna and W Bauwens, International Journal of Climatology, 2012, “A stochastic space-time model for the generation of daily rainfall in the Gaza Strip”
          # fix the correl matrix
          newMat <- a
          cholError <- TRUE
          iter <- 0
          while (cholError) {

            iter <- iter + 1
            # cat("iteration ", iter, "\n")

            # replace -ve eigen values with small +ve number
            newEig <- eigen(newMat)
            newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)

            # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
            # eig vectors
            newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

            # normalize modified matrix eqn 6 from Brissette et al 2007
            newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))

            # try chol again
            cholStatus <- try(u <- chol(newMat), silent = TRUE)
            cholError <- ifelse(inherits(cholStatus, "try-error"), TRUE, FALSE)
          }
          root <- cholStatus

          colnames(root) <- rownames(root) <- colnames(mathessian)

          pseudoV <- root %*% t(root)

          res <- sqrt(diag(pseudoV))
        }
      }
    }
  }


  SEInf <- namesp[!namesp %in% names(res)]

  res <- c(res, structure(rep(+Inf, length(SEInf)), .Names=SEInf))

  if (hessian) {
    return(list(SE=res, hessian=mathessian))
  } else {
    return(res)
  }

}
################################################################################


################################################################################
#### script1_funcoes_ULARMA


##densidade da ULARMA
dUL <- function(y,mu){
  #padronização em termos da média
  d <- ((1-mu)^2/(mu*(1-y)^3))* exp(-((y*(1-mu))/(mu*(1-y))))
  return(d)
}

##distribuição acumulada
pUL <- function(y,mu){
  p <- 1 - ((mu*y-1)/(y-1)) * exp(-((y*(1-mu))/(mu*(1-y))))
  return(p)
}

##função quantil
qUL <- function(y,mu){
  q <- ((1/mu) + pracma::lambertWn((1/mu)*(y-1)*exp(-(1/mu))))/(1 + pracma::lambertWn((1/mu)*(y-1)*exp(-(1/mu))))
  return(q)
}

##gerando valores da distribuição
rUL=function(n,mi){
  a=1/mi-1
  require(LindleyR)
  x=rlindley(n, a, mixture=TRUE)
  z=x/(1+x)
  return(z)
}


 ##############################################################################################
 #### script3_ajuste_ULARMA

#setwd("C:/Users/santa/OneDrive/Área de Trabalho/Mestrado/Dissertação/Scripts")

ULARMA.fit<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid){
  #y: série modelada
  #ar: número correspondente ao termo autorregressivo
  #ma: número corresposdnente ao termo de médias móveis
  #link: função de ligação utilizada
  #names_phi: ordem autorregressiva
  #names_theta: ordem médias móveis
  #names_beta: nome das variáveis exógenas
  #diag: gráficos (0 - sem gráficos, 1 - gráficos de diagnósticos, 2 - gráficos de diagnósticos e pdf)
  #h1: número de previsões realizadas
  #X: matrix de covariáveis exógenas
  #X_hat: matriz de covariáveis exógenas para previsões
  #resid: define o tipo de resíduo utilizado (1 - padronizado, 2 - deviance, 3 - quartílico)

  #pacotes necessários

  #definindo o número máximo de iterações para atingir a convergência
  maxit1<-100

  #saida/resultado da função
  z <- c()

  #ajeitando as saídas
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink

  #aplicando a função de ligação na resposta
  ynew = linkfun(y)
  ystar = log(y/(1-y))#faz pra comparação com codigo acima??

  #Identificando os termos temporais
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)

  #vetor de previsões
  y_prev <- c(rep(NA,(n+h1)))

  # Estimando as covariáveis do termo autorregressivo
  if(any(is.na(ar)==F)){
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)#matriz com a depêndencia temporal de ordem p

    for(i in 1:(n-m)){
      P[i,] <- ynew[i+m-ar]
    }

    Z <- cbind(rep(1,(n-m)),P)
  } else{
    Z <- as.matrix(rep(1,(n-m)))
  }

  # matriz de desenho do modelo (com ou sem variáveis exógenas)
  if(any(is.na(X)==T)){# sem variávevis exógenas
    x <- as.matrix(Z)
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    ajuste = lm.fit(x, Ynew)#ajuste do modelo linear para estimar a dependência
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)

  }else{
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n1 - k) * (dlink)^2)# pq aparece o termo dlink no denominador? (aula 6 de MLG)
  }

  ######### Estrutura do modelo sem a presença de covariáveis exógenas, somente a estrutura temporal (ARMA)
  ##neste caso é considerado os termos autorregressivos e de médas móveis
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T)){
    reg <- c(mqo, rep(0,q1)) # valores iniciais para os parâmetros

    #função de log-verossimilhança
    loglik <- function(z){
      #z é o vetor de parâmetros
      #alpha é constante
      #phi é termos autorregressivos
      #theta é medias móveis
      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]

      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_phi,names_theta)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }


    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+1)]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]

    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
  }

  ##neste caso é considerado somente o termo autorregressivo
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T)){

    q1<-0
    reg <- c(mqo) # valores iniciais dos parâmetros

    #função de log-verossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]

      eta<-rep(NA,n)

      for(i in (m+1):n){
        eta[i]<-alpha + (phi%*%ynew[i-ar])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_phi)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }

    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+1)]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]

    z$alpha <- alpha
    z$phi <- phi

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + (phi%*%ynew[i-ar])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
    }

  }

  ##neste caso é considerado somente o termo médias móveis
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T)){
    p1<-0
    reg <- c(mqo,rep(0,q1)) # valores iniciais dos parâmetros

    loglik <- function(z){
      alpha <- z[1]
      theta <- z[2:(q1+1)]

      eta <- error <- rep(0,n) # E(error)=0

      for(i in (m+1):n){
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<-ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_theta)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }

    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+1)]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    theta <- coef[2:(q1+1)]

    z$alpha <- alpha
    z$theta <- theta

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + (theta%*%errorhat[i-ma])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }

  ######### Estrutura do modelo com a presença de covariáveis exógenas
  ##neste caso é considerado os termos temporais no modelo (com covariáveis)
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F)){
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1),beta1) # valores iniciais dos parâmetros

    #função de logverossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]
      beta <- z[(p1+q1+2):length(z)]

      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)

      for(i in (m+1):n){

        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] # predictor scale
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_phi,names_theta,names_beta)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }

    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+1+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    beta <- coef[(p1+q1+2):length(coef)]


    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)


    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
  }

  ##neste caso é considerado somente os termos autorregressivos no modelo (com covariáveis)
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F)){
    q1<-0
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], beta1) # valores iniciais dos parâmetros

    #função de logverossimilhança
    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]
      beta <- z[(p1+2):length(z)]

      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_phi,names_beta)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }

    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+1+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    beta <- coef[(p1+2):length(coef)]

    z$alpha <- alpha
    z$phi <- phi

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)


    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
  }

  ##neste caso é considerado somente os termos de medias moveis no modelo (com covariáveis)
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F)){
    p1<-0
    beta1<- mqo[(2):length(mqo)]
    reg <- c(mqo[1], rep(0,q1), beta1) # valores iniciais dos parâmetros

    #função de log-verossimilhança
    loglik <- function(z){
      alpha <- z[1]
      theta = z[(2):(q1+1)]
      beta <- z[(q1+2):length(z)]

      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]

      ll <- suppressWarnings( log( dUL(y1,mu) ) )
      sum(ll)
    }

    names_par <- c("alpha",names_theta,names_beta)

    opt <- optim(reg, loglik,
                 method = "BFGS",
                 control = list(fnscale = -1),
                 hessian = TRUE)

    if (opt$conv != 0){
      warning("THE FUNCTION DID NOT CONVERGE!")
    }
    #else{
    #   warning("A FUNÇÃO CONVERGIU!")
    # }

    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+1+ncol(X) )]
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    theta <- coef[(2):(q1+1)]
    beta <- coef[(q1+2):length(coef)]

    z$alpha <- alpha
    z$theta <- theta

    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    K <- opt$hessian
    z$K <- K

    #### Previsões
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)#X[i - ar, ]


    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }

  }

  z$serie <- y
  z$ularma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]

  ########################################## residuals
  res1 <- (y-z$fitted) #resíduo ordinário

  vary <- z$fitted[(m+1):n]*(((1/z$fitted[(m+1):n])-1)^2*exp((1/z$fitted[(m+1):n])-1) *expint_En(x = ((1/z$fitted[(m+1):n])-1), order = 1)
                             -(1/z$fitted[(m+1):n])+2) - z$fitted[(m+1):n]^2

  #resíduo pradronizado
  z$resid1 <- (res1[(m+1):n]/sqrt(vary))

  #cálculo da deviance
  l_tilde <- log(dUL(y,y))
  l_hat <- log(dUL(y,z$fitted))

  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)]<-0

  z$l_hat <- l_hat

  #resíduo da deviance
  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))

  #resíduos quartilico
  z$resid3 <- as.vector(qnorm(pUL(y[(m+1):n],z$fitted[(m+1):n])))


  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3

  #matriz de informação observada
  #aqui é feita a decomposição de cholesk para verificar se a matriz é positiva definida, ou seja, inversível (https://www.ime.unicamp.br/~marcia/AlgebraLinear/Arquivos%20PDF/demo_cholesky.pdf)
  Kchol<- tryCatch(chol(z$K), error = function(e) return("error"),
                   warning = function(o) return("error"))#https://r-lang.com/r-trycatch-function/

  if(Kchol[1] == "error"){
    z$vcov <- try(solve(K))#https://blog.curso-r.com/posts/2017-04-09-try/
    warning("We have problems with information matrix inversion!")

  }else{
    vcov <- try(chol2inv(Kchol))
    z$vcov <- vcov
  }

  #desvio padrão
  stderror <- sqrt(diag(z$vcov))
  if(any(is.na(stderror))){
    stderror <- tryCatch(SEfromHessian(opt$hessian),
                         error = function(e){
                           #library(numDeriv)

                           rep(1,length(z$coef))
                         } ,
                         warning = function(o) return("error"))

  }
  z$stderror <- stderror

  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )

  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])

  #medidas comparação
  if(any(is.na(X)==F)){
    z$k<- (p1+q1+1+length(beta))
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+1)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }

  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model <- model_presentation
  z$link <- link

  # if(diag>0){
  #
  #   print(model_presentation)
  #   print(" ",quote=F)
  #   print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
  #   print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
  #   print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
  #
  #   print("Residuals:",quote=F)
  #   print(summary(residc))
  #
  #   t<-seq(-5,n+6,by=1)
  #
  #   par(mfrow=c(1,1))
  #   par(mar=c(2.8, 2.7, 1.2, 1))
  #   par(mgp=c(1.7, 0.45, 0))
  #   plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
  #   lines(t,rep(-3,n+12),lty=2,col=1)
  #   lines(t,rep(3,n+12),lty=2,col=1)
  #   lines(t,rep(-2,n+12),lty=3,col=1)
  #   lines(t,rep(2,n+12),lty=3,col=1)
  #
  #
  #   max_y<- max(c(z$fitted,y),na.rm=T)
  #   min_y<- min(c(z$fitted,y),na.rm=T)
  #   plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
  #        xlab="Fitted values",ylab="Observed data",
  #        xlim=c(0.95*min_y,max_y*1.05),
  #        ylim=c(0.95*min_y,max_y*1.05))
  #   lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
  #
  #   plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
  #        xlab="Fitted values",ylab="Residuals")
  #
  #   densidade<-density(residc)
  #   plot(densidade,ylab="density",main=" ")
  #   lines(densidade$x,dnorm(densidade$x),lty=2)
  #   legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
  #          pt.bg="white", lty=c(1,2), bty="n")
  #
  #   acf(residc,ylab="ACF",xlab="Lag")
  #
  #   pacf(residc,ylab="PACF",xlab="Lag")
  #
  #   max_r<- max(residc,na.rm=T)
  #   min_r<- min(residc,na.rm=T)
  #   qqnorm(residc, pch = "+",
  #          xlim=c(0.95*min_r,max_r*1.05),
  #          ylim=c(0.95*min_r,max_r*1.05),
  #          main="",xlab="Normal quantiles",ylab="Empirical quantiles")
  #   lines(c(-10,10),c(-10,10),lty=2)
  #
  #   par(mfrow=c(1,1))
  #   plot(y,type="l",ylab="serie",xlab="tempo")
  #   lines(z$fitted,col="red")
  #
  #   fim<-end(y)[1]+end(y)[2]/12
  #
  #   y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
  #   par(mfrow=c(1,1))
  #   plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
  #   abline(v=fim,lty=2)
  #   lines(y)
  #
  #   w1<-5
  #   h1<-4
  #
  #   if(diag>1)
  #   {
  #     png(file = "resid_v_ind.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
  #       lines(t,rep(-3,n+12),lty=2,col=1)
  #       lines(t,rep(3,n+12),lty=2,col=1)
  #       lines(t,rep(-2,n+12),lty=3,col=1)
  #       lines(t,rep(2,n+12),lty=3,col=1)
  #     }
  #     dev.off()
  #
  #     png(file = "resid_v_fitted.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
  #            xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
  #       lines(t,rep(-3,n+12),lty=2,col=1)
  #       lines(t,rep(3,n+12),lty=2,col=1)
  #       lines(t,rep(-2,n+12),lty=3,col=1)
  #       lines(t,rep(2,n+12),lty=3,col=1)
  #     }
  #     dev.off()
  #
  #     png(file = "obs_v_fit.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
  #            xlab="Fitted values",ylab="Observed data",
  #            xlim=c(0.95*min_y,max_y*1.05),
  #            ylim=c(0.95*min_y,max_y*1.05))
  #       lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
  #     }
  #     dev.off()
  #
  #     png(file = "resid_density.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(1.5, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #
  #       plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
  #       lines(densidade$x,dnorm(densidade$x),lty=2)
  #       legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
  #              pt.bg="white", lty=c(1,2), bty="n")
  #     }
  #     dev.off()
  #
  #     png(file = "resid_FAC.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       acf(residc,ylab="ACF",xlab="Lag")
  #     }
  #     dev.off()
  #
  #     png(file = "resid_FACP.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       pacf(residc,ylab="PACF",xlab="Lag")
  #     }
  #     dev.off()
  #
  #     png(file = "qq_plot.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       qqnorm(residc, pch = "+",
  #              xlim=c(0.95*min_r,max_r*1.05),
  #              ylim=c(0.95*min_r,max_r*1.05),
  #              main="",xlab="Normal quantiles",ylab="Empirical quantiles")
  #       lines(c(-10,10),c(-10,10),lty=2)
  #     }
  #     dev.off()
  #
  #     png(file = "adjusted.png",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(y,type="l",ylab="Serie",xlab="Time")
  #       lines(z$fitted,col="red")
  #     }
  #     dev.off()
  #
  #
  #     png(file = "forecast.png",horizontal=F,paper="special",width = 6, height = 4.7,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(y_prev,type="l",lty=2,col="red", ylim=c(min(y),max(y)),ylab="RH",xlab="Times")
  #       abline(v=fim,lty=2)
  #       lines(y)
  #       legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
  #              pt.bg="white", lty=c(1,2), bty="n",col=c(1,"red"))
  #     }
  #     dev.off()
  #
  #   }
  # }

  return(z)



}
#################################################################################






################################################################################
#### script4_geracao_valores_v2

#setwd("C:/Users/santa/OneDrive/Área de Trabalho/Mestrado/Dissertação/Scripts")


simu.ularma_v2 <- function(n,phi=NA,theta=NA,alpha=0.0,link="logit",X=NA,beta=NA)
{
  ar<-NA
  ma<-NA

  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }

  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }

  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  }

  link <- structure(list(link = linktemp,
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )

  linkfun <- link$linkfun
  linkinv <- link$linkinv

  if(any(is.na(X) == TRUE)){
    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F))
    {
      #print("ARMA model")
      # seasonal part

      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)

      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }

      return(ts(y[(m+1):(n+m)]))

    } # ARMA model


    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T))
    {
      #print("AR model")

      p <- max(ar)
      m <- 2*p

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      eta<- y <- rep(NA,n+m)

      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + (phi%*%ynew[i-ar])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
      }

      return( ts(y[(m+1):(n+m)]))
    } # AR model


    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F))
    {
      #print("MA model")

      q <- max(ma)
      m <- 2*q

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      eta <- y <- error <- rep(0,n+m) # E(error)=0

      #print(ma)

      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }

      return( ts(y[(m+1):(n+m)]))
    } # fim MA model

  } else {

    ###### ARMA model
    if(any(is.na(phi)==F) && any(is.na(theta)==F))
    {
      #print("ARMA model")
      # seasonal part

      p <- max(ar)
      q <- max(ma)
      m <- 2*max(p,q)

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      error<-rep(0,n+m) # E(error)=0
      eta<- y <- rep(NA,n+m)

      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)

      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + as.numeric(X[i,]%*%as.matrix(beta)) + as.numeric((phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))) + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]

      }

      return(ts(y[(m+1):(n+m)]))

    } # ARMA model


    ###### AR model
    if(any(is.na(phi)==F) && any(is.na(theta)==T))
    {
      #print("AR model")

      p <- max(ar)
      m <- 2*p

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      eta<- y <- rep(NA,n+m)
      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)
      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
      }

      return( ts(y[(m+1):(n+m)]))
    } # AR model


    ###### MA model
    if(any(is.na(phi)==T) && any(is.na(theta)==F))
    {
      #print("MA model")

      q <- max(ma)
      m <- 2*q

      ynew <-rep(alpha,(n+m))
      mu <- linkinv(ynew)

      eta <- y <- error <- rep(0,n+m) # E(error)=0
      X <- rbind(matrix(X[1:m,],ncol=length(beta)),X)
      #print(ma)

      for(i in (m+1):(n+m))
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        mu[i] <- linkinv(eta[i])
        y[i] <- rUL(1,mu[i])
        ynew[i] <- linkfun(y[i])
        error[i]<- ynew[i]-eta[i]
      }

      return( ts(y[(m+1):(n+m)]))
    } # fim MA model
  }


}

################################################################################



# @param diag
# If \code{diag = 0}: without graphs (useful for simulations)
# If \code{diag = 1}: plot diagnostic graphs
# If \code{diag = 2}: make pdf with diagnostic graphs








#' ULARMA Model
#'
#' Fit and forecast the Unit-Lindley Autoregressive and Moving Average (ULARMA) model as
#' an extension of the unit-Lindley regression for the case of autocorrelated residuals.
#' @param y   a numeric vetor or time series with UnitLindley distribution
#' @param ar    autoregressive order
#' @param ma    moving avarege order
#' @param link  the link function ("logit", "probit", "cloglog")
#' @param h forecast horizon
#' @param X optionally numerical vector or matrix of external regressors, which must have the same number of rows as y.
#' @param X_hat the future values of \code{X} to forecast, which must have \code{h} rows.
#' @param resid type of residuals (1 - standardized, 2 - deviance, 3 - quartile)
#' @returns return a list with several informations about the model
#' @examples
#' x = rbeta(50, 2, 2)
#' fit = ULARMA(ts(x), 1, 1)
#' fit$fitted
#' @author José G S Sena, Fabio M Bayer, Paulo H. Ferreira, José A Fiorucci
#' @import "expint","numDeriv"
#' @export
ULARMA <- function (y, ar=NA, ma=NA,link="logit",h=6,X=NA,X_hat=NA,resid=3){

  # "expint","dplyr","openxlsx","lubridate","forecast","ggplot2","gridExtra","ggpubr","tidyr","corrplot","numDeriv"

  #verificando se a variável possui valores dentro do intervalo unitário padrão
  if (min(y) <= 0 || max(y) >= 1){
    stop("Out of range (0,1)!")
  }

  #verificando se a variável é uma série temporal
  if(is.ts(y)==T){
    freq<-frequency(y)
  }else stop("y must be a time-series object")

  #verificando se o vetor com os termos AR/MA estão completos
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))

  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))

  #verificando a presença de covariaveis
  if(any(is.na(X))==F) names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))

  #pegando as ordens dos termos temporais
  p <- max(ar)
  q <- max(ma)
  n <- length(y)

  m <- max(p,q,na.rm=T)

  p1 <- length(ar)
  q1 <- length(ma)

  #ajeitando as saidas da função de ligação
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)#converte expressão em character
    if (linktemp == "link")
      linktemp <- eval(link)#verifica a expressão
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link function not available, use \"logit\", ",
                  "\"probit\" and \"cloglog\""))

  link1 <- structure(list(link = linktemp,
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv,
                          mu.eta = stats$mu.eta,
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t))) # ver esse ponto aqui
  )
  )


  fit1 <- ULARMA.fit(y=y, ar=ar, ma=ma, link=link1, names_phi=names_phi,
                     names_theta=names_theta, names_beta=names_beta,
                     diag=0, h1=h, X=X, X_hat=X_hat,resid=resid) # model estimation

  return(fit1)
}
###################################################################################


