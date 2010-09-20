flowClust<-function(x, expName="Flow Experiment", varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1, nu.est=0, trans=1, min.count=10, max.count=10, min=NULL, max=NULL, level=0.9, u.cutoff=NULL, z.cutoff=0, randomStart=10, B.init=B, tol.init=1e-2, seed=1, criterion="BIC", control=NULL) 
{
    if (is(x, "flowFrame")) {
        if (length(varNames)==0) {
            y <- exprs(x)
            varNames <- colnames(y)
        }
        else {
            y <- as.matrix(exprs(x)[, varNames])
        }
    }
    else if (is(x, "matrix")) {
        if (length(varNames)==0) {
            y <- x
            if (length(colnames(x))==0) 
                varNames <- "Not Available"
            else varNames <- colnames(x)
        }
        else {
            y <- as.matrix(x[, varNames])
        }
    }
    else if (is(x, "data.frame")) {
        if (length(varNames)==0) {
            y <- as.matrix(x)
            varNames <- colnames(x)
        }
        else {
            y <- as.matrix(x[, varNames])
        }
    }
    else if (is(x, "vector")) {
        y <- matrix(x)
        if (length(varNames)==0) 
            varNames <- "Not Available"
    }
    else {
        stop(paste("Object ", as.character(x), " is not of class flowFrame / matrix / data frame!"))
    }


    # finding filtered observations
    rm.max <- rm.min <- rep(FALSE, nrow(y))
    if (max.count > -1) {
        if (is.null(max)[1]) 
            max <- apply(y, 2, max)
        for (k in 1:ncol(y))  if (sum(y[,k]>=max[k]) >= max.count) 
            rm.max <- rm.max | (y[,k] >= max[k])
    }
    if (min.count > -1) {
        if (is.null(min)[1]) 
            min <- apply(y, 2, min)
        for (k in 1:ncol(y))  if (sum(y[,k]<=min[k]) >= min.count) 
            rm.min <- rm.min | (y[,k] <= min[k])
    }
    include <- !rm.max & !rm.min
    
    if((length(grep("multicore",loadedNamespaces()))==0) & (length(grep("snowfall",loadedNamespaces()))==0 || suppressMessages(!sfParallel())))
    {
      message("Using the serial version of flowClust")    
      # C version
      result<-lapply(as.list(1:length(K)),.flowClustK, y, expName=expName, varNames=varNames, K=K, B=B, tol=tol, nu=nu, lambda=lambda, nu.est=nu.est, trans=trans, min.count=min.count, max.count=max.count, min=min, max=max, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff, randomStart=randomStart, B.init=B.init, tol.init=tol.init, seed=seed, criterion=criterion, control=control,include=include, rm.max, rm.min)
    }
    else if(length(grep("multicore",loadedNamespaces()))==1)
    {
      cores<-getOption("cores")
      if(is.null(cores))
      {
        nClust<-multicore:::volatile$detectedCores
      }
      else
      {
        nClust<-cores
      }
      message("Using the multicore version of flowClust with ",nClust," cores")
      # Split into nClust segReadsList
      result<-mclapply(as.list(1:length(K)),.flowClustK, y, expName=expName, varNames=varNames, K=K, B=B, tol=tol, nu=nu, lambda=lambda, nu.est=nu.est, trans=trans, min.count=min.count, max.count=max.count, min=min, max=max, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff, randomStart=randomStart, B.init=B.init, tol.init=tol.init, seed=seed, criterion=criterion, control=control,include=include, rm.max, rm.min, mc.preschedule=FALSE)
    }
    else if(length(grep("snowfall",loadedNamespaces()))==1 && sfParallel())
    {
      # Number of clusters
      nClust<-sfCpus()
      message("Using the parallel (snowfall) version of flowClust with ", nClust, " cpus or cores")
      result<-sfLapply(as.list(1:length(K)),.flowClustK, y, expName=expName, varNames=varNames, K=K, B=B, tol=tol, nu=nu, lambda=lambda, nu.est=nu.est, trans=trans, min.count=min.count, max.count=max.count, min=min, max=max, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff, randomStart=randomStart, B.init=B.init, tol.init=tol.init, seed=seed, criterion=criterion, control=control,include=include, rm.max, rm.min)
    }

    # Simply return a flowClust object
    if (length(K)==1)
    {
        result[[1]]
    }
    # Create a list flowClustList
    else
    {
        result <- new("flowClustList", result, criterion=criterion)
        result@index <- which.max(criterion(result, criterion))
        result
    }
}

.flowClustK<-function(i, y, expName="Flow Experiment", varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1, nu.est=0, trans=1, min.count=10, max.count=10, min=NULL, max=NULL, level=0.9, u.cutoff=NULL, z.cutoff=0, randomStart=10, B.init=B, tol.init=1e-2, seed=1, criterion="BIC", control=NULL, include, rm.max,rm.min)
{
  y <- as.matrix(y[include,])
  ly <- nrow(y)
  py <- ncol(y)
  if (min(y)<=0 && lambda<=0) 
      stop("lambda must be positive when data contain zero / negative values!")


  # to determine the rule of calling outliers
  if (nu != Inf) {
      if (is.null(u.cutoff)) {
          if (!nu.est) {
              cc <- py * qf(level, py, nu)
              u.cutoff <- (nu + py) / (nu + cc)
          }
          else {
              u.cutoff <- level     # pass level to .C call to compute u.cutoff
          }
          ruleOutliers <- c(0, level, z.cutoff)     # 0 means quantile
      }
      else {
          ruleOutliers <- c(1, u.cutoff, z.cutoff)     # 1 means cutoff
      }
  }
  else {
      if (level != 1)
          q.cutoff <- qchisq(level, py)
      else
          q.cutoff <- -1     # -1 means no outlier identification
      ruleOutliers <- c(0, level, z.cutoff)
  }
  

  if (is.null(control$B.lambda)) control$B.lambda <- B    # BSolve=100
  if (is.null(control$B.brent)) control$B.brent <- 10000    # iterSolveMax=50
  if (is.null(control$tol.brent)) control$tol.brent <- 1e-5    # DiffSolve=1e-3
  if (is.null(control$xLow)) control$xLow <- 0.1    # xLow=.1
  if (is.null(control$xUp)) control$xUp <- 10    # xUp=1
  if (is.null(control$nuLow)) control$nuLow <- 2    # nuLow=2
  if (is.null(control$nuUp)) control$nuUp <- 100    # nuUp=30
  if (is.null(control$seed)) control$seed <- TRUE
  
  ind <- 0
  if (K[i]==1)
  { 
    label <- rep(1, ly)
  }
  else if (!randomStart)
  {
    if (py==1)
    {
      q <- quantile(y, seq(from=0, to=1, by=1/K[i]))
      label <- rep(0, ly)
      q[1] <- q[1]-1
      for (k in 1:K[i]) label[y>q[k] & y<=q[k+1]] <- k
    }
  }
  else 
  { # Initialization based on short EMs with random partitions if randomStart=TRUE
    if (control$seed) set.seed(seed)
    if (randomStart==1) 
    {
      label <- sample(1:K[i], ly, replace=T)
    }
    else 
    {
      maxLabel <- vector("list",randomStart)
      maxLogLike <- rep(NA,randomStart)
      for (j in 1:randomStart) 
      {
        label <- sample(1:K[i], ly, replace=TRUE)
        if (nu != Inf) 
        {
          obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly),
          as.integer(py), as.integer(K[i]),
          w=rep(0,K[i]), mu=rep(0,K[i]*py),
          precision=rep(0,K[i]*py*py), 
          lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
          nu=as.double(rep(nu,K[i])),
          z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
          as.integer(label), uncertainty=double(ly), 
          as.double(rep(u.cutoff,K[i])), as.double(z.cutoff), 
          flagOutliers=integer(ly), as.integer(B.init), 
          as.double(tol.init), as.integer(trans),
          as.integer(nu.est), logLike=as.double(0), 
          as.integer(control$B.lambda), as.integer(control$B.brent), 
          as.double(control$tol.brent), as.double(control$xLow), 
          as.double(control$xUp), as.double(control$nuLow), 
          as.double(control$nuUp), PACKAGE="flowClust"))
        }
        else 
        {
          obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly), 
          as.integer(py), as.integer(K[i]), 
          w=rep(0,K[i]), mu=rep(0,K[i]*py), 
          precision=rep(0,K[i]*py*py),
          lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))), 
          z=rep(0, ly*K[i]), u=rep(0,ly*K[i]),
          as.integer(label), uncertainty=double(ly), 
          as.double(q.cutoff), as.double(z.cutoff), 
          flagOutliers=integer(ly), as.integer(B.init), 
          as.double(tol.init), as.integer(trans), 
          logLike=as.double(0),
          as.integer(control$B.lambda), as.integer(control$B.brent), 
          as.double(control$tol.brent), as.double(control$xLow), 
          as.double(control$xUp),
          PACKAGE="flowClust"))
        }
        if (class(obj)!="try-error") 
        {
          maxLabel[[j]] <- label
          maxLogLike[j] <- obj$logLike
        }
      }
      ind <- order(maxLogLike, decreasing=T, na.last=NA)
    }
  }

  # long EMs
  for (M in ind) 
  {
    if (nu != Inf) 
    {
      obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly), 
      as.integer(py), as.integer(K[i]),
      w=rep(0,K[i]), mu=rep(0,K[i]*py),
      precision=rep(0,K[i]*py*py),
      lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))), 
      nu=as.double(rep(nu,K[i])),
      z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
      as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly),
      as.double(rep(u.cutoff,K[i])), as.double(z.cutoff),
      flagOutliers=integer(ly), as.integer(B),
      as.double(tol), as.integer(trans), 
      as.integer(nu.est), logLike=as.double(0),
      as.integer(control$B.lambda), as.integer(control$B.brent), 
      as.double(control$tol.brent), as.double(control$xLow), 
      as.double(control$xUp), as.double(control$nuLow),
      as.double(control$nuUp), PACKAGE="flowClust"))
    }
    else 
    {
      obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly), 
      as.integer(py), as.integer(K[i]),
      w=rep(0,K[i]), mu=rep(0,K[i]*py),
      precision=rep(0,K[i]*py*py),
      lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))), 
      z=rep(0,ly*K[i]), u=rep(0,ly*K[i]), 
      as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly), 
      as.double(q.cutoff), as.double(z.cutoff),
      flagOutliers=integer(ly), as.integer(B),
      as.double(tol), as.integer(trans), 
      logLike=as.double(0),
      as.integer(control$B.lambda), as.integer(control$B.brent), 
      as.double(control$tol.brent), as.double(control$xLow), 
      as.double(control$xUp),
      PACKAGE="flowClust"))
      if (class(obj)!="try-error")
      obj$nu <- Inf
    }
    if (class(obj)!="try-error")
    break
  }
  if (class(obj)=="try-error") 
  stop(geterrmessage())

  # output obj$precision to sigma
  sigma <- array(0, c(K[i], py, py))
  precision <- matrix(obj$precision, K[i], py * py, byrow=TRUE)
  for (k in 1:K[i])
  sigma[k,,] <- matrix(precision[k,], py, py, byrow = TRUE)

  # output BIC & ICL
  BIC <- 2*obj$logLike - log(ly) * (K[i]*(py+1)*py/2 + K[i]*py + K[i]-1 + (if (trans>1) K[i] else trans) + (if (nu.est>1) K[i] else abs(nu.est)))
  z <- matrix(obj$z, ly, K[i], byrow = TRUE)
  ICL <- BIC + 2 * sum(z*log(z), na.rm = TRUE)

  # output z, u, label, uncertainty, flagOutliers
  z <- u <- matrix(NA, length(include), K[i])
  z[include,] <- matrix(obj$z, ly, K[i], byrow=TRUE)
  u[include,] <- matrix(obj$u, ly, K[i], byrow=TRUE)
  tempLabel <- if (M==0) label else maxLabel[[M]]
  label <- uncertainty <- flagOutliers <- rep(NA, length(include))
  label[include] <- tempLabel
  uncertainty[include] <- obj$uncertainty
  flagOutliers[include] <- as.logical(obj$flagOutliers)

  result<- new("flowClust", expName=expName, varNames=varNames, K=K[i],
  w=obj$w, mu=matrix(obj$mu, K[i], py, byrow=TRUE), sigma=sigma,
  lambda=(if (trans>0) obj$lambda else if (lambda!=1) lambda else numeric(0)), nu=(if (nu.est>1) obj$nu else obj$nu[1]), z=z,
  u=u, label=label, uncertainty=uncertainty, 
  ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min), 
  rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL)
  result
}
