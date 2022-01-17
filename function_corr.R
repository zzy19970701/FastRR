cor.part <- function(Xdata,Ydata,pvc=0.01){
  
  ## regist Multicore
  library(doParallel)
  cl.cores <- detectCores()
  if (cl.cores <= 2) {
    cl.cores <- 1}
  if (cl.cores > 2) {
    if (cl.cores > 10) {
      cl.cores <- 10
    }
    else {
      # cl.cores <- detectCores() - 1
      cl.cores <- 5
    }
  }
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  ###########
  # Xdata <- G_c
  # Ydata <- Y_c
  
  le <- ncol(Xdata)
  corv <- vector()
  corp <- vector()
  ##correlation
  corallres = foreach(i = 1:le, .multicombine = TRUE, .combine = "rbind") %dopar% {
    if (var(Xdata[, i]) > 0) {
      corp[i] <- cor.test(Xdata[, i], Ydata)$p.value
      corv[i] <- abs(cor(Xdata[, i], Ydata))
      corall <- c(corp[i], corv[i])
    }
    else {
      corp[i] <- 1
      corv[i] <- 0
      corall <- c(corp[i], corv[i])
    }
  }
  rownames(corallres) <- NULL
  corp <- corallres[, 1]
  corv <- corallres[, 2]
  stopCluster(cl)
  rm(Xdata)
  gc()
  
  if (length(which(corp < pvc)) <= nrow(y)) {
    ee <- as.vector(which(corp < pvc))
  }
  else {
    n1 <- nrow(Ydata) - 1
    ee <- as.vector(which(rank(corv) >= (le - n1), arr.ind = T))
  }
  return(ee)
}






