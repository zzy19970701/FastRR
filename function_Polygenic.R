{
  emma.eigen.L <- function(Z,K,complete=TRUE) {
    if ( is.null(Z) ) {
      return(emma.eigen.L.wo.Z(K))
    }
    else {
      return(emma.eigen.L.w.Z(Z,K,complete))
    }
  }
  
  emma.eigen.L.wo.Z <- function(K) {
    eig <- eigen(K,symmetric=TRUE)
    return(list(values=eig$values,vectors=eig$vectors))
  }
  
  emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
    return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
  }
  
  emma.eigen.R <- function(Z,K,X,complete=TRUE) {
    if ( ncol(X) == 0 ) {
      return(emma.eigen.L(Z,K))
    }
    else if ( is.null(Z) ) {
      return(emma.eigen.R.wo.Z(K,X))
    }
    else {
      return(emma.eigen.R.w.Z(Z,K,X,complete))
    }
  }
  
  emma.eigen.R.wo.Z <- function(K, X) {
    n <- nrow(X)
    q <- ncol(X)
    S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
    eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
    stopifnot(!is.complex(eig$values))
    return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
  }
  
  emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
    if ( complete == FALSE ) {
      vids <-  colSums(Z) > 0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    
    SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
    eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
    if ( is.complex(eig$values) ) {
      eig$values <- Re(eig$values)
      eig$vectors <- Re(eig$vectors)    
    }
    qr.X <- qr.Q(qr(X))
    return(list(values=eig$values[1:(t-q)],
                vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                             complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
  }
  
  emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
  }
  
  emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
    
  }
  
  emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
    n <- length(xi)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
  }
  
  emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
  }
  
  emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <-  exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
  }
  
  emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <-  exp(logdelta)
    return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
  }
  
  emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
    nq <- length(etas)
    delta <- exp(logdelta)
    etasq <- etas*etas
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
  }
  
  emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
    t <- t1
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    etasq <- etas.1*etas.1
    ldelta <- delta*lambda+1
    return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
  }
  
  emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL)
  {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      warning("X is singular")
      return (list(ML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.wo.Z(K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m)    
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
      Xis.1<-matrix(eig.L$values,n,m)
      Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1 
      Etasq <- matrix(etas*etas,n-q,m)
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
    }
    else {
      if ( is.null(eig.L) ) {
        eig.L <- emma.eigen.L.w.Z(Z,K)
      }
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,t-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      
      Xis.1<-matrix(eig.L$values,t,m)
      Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1  
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
    }
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
    }
    maxvg <- maxve*maxdelta
    
    return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  }
  
  emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL) {
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    
stopifnot(ncol(K)==t)
stopifnot(nrow(X)==n)
    
    if ( det(crossprod(X,X)) == 0 ) {
      warning("X is singular")
      return (list(REML=0,delta=0,ve=0,vg=0))
    }
    
    if ( is.null(Z) ) {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.wo.Z(K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1<-matrix(eig.R$values,n-q,m)
      Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas*etas,n-q,m)
      
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
      
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
    }
    else {
      if ( is.null(eig.R) ) {
        eig.R <- emma.eigen.R.w.Z(Z,K,X)
      }
      etas <- crossprod(eig.R$vectors,y)
      etas.1 <- etas[1:(t-q)]
      etas.2 <- etas[(t-q+1):(n-q)]
      etas.2.sq <- sum(etas.2*etas.2)
      
      logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
      m <- length(logdelta)
      delta <- exp(logdelta)
      
      Lambdas.1 <- matrix(eig.R$values,t-q,m) 
      Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
      Etasq <- matrix(etas.1*etas.1,t-q,m)
      
      dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
      
      optlogdelta <- vector(length=0)
      optLL <- vector(length=0)
      if ( dLL[1] < esp ) {
        optlogdelta <- append(optlogdelta, llim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      if ( dLL[m-1] > 0-esp ) {
        optlogdelta <- append(optlogdelta, ulim)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
      }
      
      for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
    }  
    
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    
    if ( is.null(Z) ) {
      maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
    }
    else {
      maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
    }
    maxvg <- maxve*maxdelta
    return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  }
  
  emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
    if( is.null(Z) ){
      return(emma.maineffects.B.Zo(K,deltahat.g))
    }
    else{
      return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
    }
  }
  
  emma.maineffects.B.Zo <-function(K,deltahat.g){
    t <- nrow(K)
    stopifnot(ncol(K) == t)
    
    B<-deltahat.g*K+diag(1,t)
    eig<-eigen(B,symmetric=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A))
  }
  
  emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
    if ( complete == FALSE ) {
      vids <- colSums(Z)>0
      Z <- Z[,vids]
      K <- K[vids,vids]
    }
    
    n <- nrow(Z)  
    B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
    eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
    qr.B<-qr(B)
    q<-qr.B$rank
    
    stopifnot(!is.complex(eig$values))
    
    A<-diag(1/sqrt(eig$values[1:q]))
    Q<-eig$vectors[,1:q]
    C<-Q%*%A%*%t(Q)
    return(list(mC=C,Q=Q,A=A,complete=TRUE))
  }
  
  emma.MLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    
    stopifnot(nrow(W_c)==n)
    
    M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    etas<-crossprod(M_c,Y_c)
    
    LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
    return(list(ML=LL))
    
  }
  
  emma.REMLE0.c <- function(Y_c,W_c){
    
    n <- length(Y_c)
    
    stopifnot(nrow(W_c)==n)
    
    M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
    eig <-eigen(M_c)
    t <-qr(W_c)$rank
    v <-n-t
    U_R <-eig$vector[,1:v]
    etas<-crossprod(U_R,Y_c)
    
    LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
    return(list(REML=LL))
  }
  
  multinormal<-function(y,mean,sigma)
  {
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  
}