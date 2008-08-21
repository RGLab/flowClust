flowClust<-function(x, expName="Flow Experiment", varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1, trans=TRUE, min.count=10, max.count=10, min=NULL, max=NULL, level=0.9, u.cutoff=NULL, z.cutoff=0, randomStart=FALSE, seed=1)
{

    if(is(x,"flowFrame"))
    {
        if(length(varNames)==0)
        {
            y<-exprs(x)
            varNames<-colnames(y)
        }
        else
        {
            y<-as.matrix(exprs(x)[,varNames])
        }
    }
    else if(is(x,"matrix"))
    {
        if(length(varNames)==0)
        {
            y<-x
            if (length(colnames(x))==0) varNames<-"Not Available"  else varNames<-colnames(x)
        }
        else
        {
            y<-as.matrix(x[,varNames])
        }
    }
    else if(is(x,"data.frame"))
    {
        if(length(varNames)==0)
        {
            y<-as.matrix(x)
            varNames<-colnames(x)
        }
        else
        {
            y<-as.matrix(x[,varNames])
        }
    }
    else if(is(x,"vector"))
    {
        y<-matrix(x)
        if(length(varNames)==0) varNames<-"Not Available"
    }
    else
    {
      stop(paste("Object ", as.character(x)," is not of class flowFrame / matrix / data frame!"))
    }

    rm.max <- rm.min <- rep(FALSE, nrow(y))
    if (max.count > -1)
    {
        if (is.null(max)[1]) max <- apply(y,2,max)
        for (k in 1:ncol(y))  if (sum(y[,k]>=max[k]) >= max.count)  rm.max <- rm.max | (y[,k]>=max[k])
    }
    if (min.count > -1)
    {
        if (is.null(min)[1]) min <- apply(y,2,min)
        for (k in 1:ncol(y))  if (sum(y[,k]<=min[k]) >= min.count)  rm.min <- rm.min | (y[,k]<=min[k])
    }
    include <- !rm.max & !rm.min


    y <- as.matrix(y[include,])
    ly <- nrow(y)
    py <- ncol(y)
    if (min(y)<=0 && lambda<=0) stop("lambda must be positive when data contain zero / negative values!")

    if(nu!=Inf)
    {
        if(is.null(u.cutoff)) 
        {
            cc <- py * qf(level, py, nu)
            u.cutoff <- (nu+py) / (nu+cc)
            ruleOutliers <- c(0, level, z.cutoff)     # 0 means quantile
        }  
        else  
        {
            ruleOutliers <- c(1, u.cutoff, z.cutoff)     # 1 means cutoff
        }
    }
    else
    {
        if (level!=1) q.cutoff <- qchisq(level, py)  else q.cutoff <- -1     # -1 means no outlier identification
        ruleOutliers <- c(0, level, z.cutoff)
    }
    

    # Initialization based on mclust
    # If more than 1500 observations, only use 1500 at random
    if (!randomStart)
    {
        if (py==1)
        {
            q <- quantile(y, seq(from=0, to=1, by=1/K))
            label <- rep(0, ly)
            q[1] <- q[1] - 1
            for (k in 1:K) label[ y>q[k] & y<=q[k+1] ] <- k
        }
        else
        {
            if(ly>1500)
            {
                set.seed(seed)
                ySubset <- sample(1:ly,1500)
            }
            else
            {
                ySubset<-1:ly
            }

            hcPairs <- hc((if (py>1) "VVV" else "V"), (if (lambda!=0) (sign(y[ySubset,])*abs(y[ySubset,])^lambda-1)/lambda else log(y[ySubset,])))
            label <- rep(0,ly)
            label[ySubset] <- hclass(hcPairs,K)
        }
    }
    else
    {
        set.seed(seed)
        label <- sample(1:K,ly,replace=T)
    }
    
    if (nu!=Inf)
    {
        obj<-.C("flowClust",
            as.double(t(y)),
            as.integer(ly),
            as.integer(py),
            as.integer(K),
            w=rep(0.0,K),
            mu=rep(0.0,K*py),
            precision=rep(0,K*py*py),
            lambda=as.double(lambda),
            nu=as.double(nu),
            z=rep(0.0,ly*K),
            u=rep(0.0,ly*K),
            as.integer(label),	
            uncertainty=double(ly),
            as.double(u.cutoff),
            as.double(z.cutoff),
            flagOutliers=integer(ly),
            as.integer(B),
            as.double(tol),
            as.integer(trans),
            logLike=as.double(0),
            package="flowClust")
    } else
    {
        obj<-.C("flowClustGaussian",
            as.double(t(y)),
            as.integer(ly),
            as.integer(py),
            as.integer(K),
            w=rep(0.0,K),
            mu=rep(0.0,K*py),
            precision=rep(0,K*py*py),
            lambda=as.double(lambda),
            z=rep(0.0,ly*K),
            u=rep(0.0,ly*K),
            as.integer(label),	
            uncertainty=double(ly),
            as.double(q.cutoff),
            as.double(z.cutoff),
            flagOutliers=integer(ly),
            as.integer(B),
            as.double(tol),
            as.integer(trans),
            logLike=as.double(0),
            package="flowClust")
    }

    # output obj$precision to sigma
    sigma<-array(0,c(K,py,py))
    precision<-matrix(obj$precision,K,py*py,byrow=TRUE)
    for(k in 1:K) sigma[k,,]<-matrix(precision[k,],py,py,byrow=TRUE)

    BIC <- 2*obj$logLike - log(ly)*(K*(py+1)*py/2+K*py+K-1+as.integer(trans))
    z<-matrix(obj$z,ly,K,byrow=TRUE)
    ICL <- BIC + 2*sum(z*log(z),na.rm=TRUE)

    z <- u <- matrix(NA, length(include), K)
    z[include,] <- matrix(obj$z,ly,K,byrow=TRUE)
    u[include,] <- matrix(obj$u,ly,K,byrow=TRUE)
    tempLabel <- label
    label <- uncertainty <- flagOutliers <- rep(NA, length(include))
    label[include] <- tempLabel
    uncertainty[include] <- obj$uncertainty
    flagOutliers[include] <- as.logical(obj$flagOutliers)

    new("flowClust", expName=expName, varNames=varNames, K=K,
        w=obj$w, mu=matrix(obj$mu,K,py,byrow=TRUE), sigma=sigma, lambda=obj$lambda,
        nu=nu, z=z, u=u, label=label, uncertainty=uncertainty,
        ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min), 
        rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL)
}
