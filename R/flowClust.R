flowClust<-function(x, expName="Flow Experiment", varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1, trans=TRUE, min.count=10, max.count=10, min=NULL, max=NULL, level=0.9, u.cutoff=NULL, z.cutoff=0, randomStart=FALSE, B.init=B, tol.init=1e-2, seed=1, criterion="BIC")
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
    

    result <- vector("list", length(K))

    for (i in 1:length(K))
    {

        if (K[i]==1) label <- rep(1,ly)
        else if (!randomStart)
        {
            if (py==1)
            {
                q <- quantile(y, seq(from=0, to=1, by=1/K[i]))
                label <- rep(0, ly)
                q[1] <- q[1] - 1
                for (k in 1:K[i]) label[ y>q[k] & y<=q[k+1] ] <- k
            }
            else
            {
                # Initialization based on mclust
                # If more than 1500 observations, only use 1500 at random
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
                label[ySubset] <- hclass(hcPairs,K[i])
            }
        }
        else
        {
            set.seed(seed)
            if (randomStart==1)
            {
                label <- sample(1:K[i],ly,replace=T)
            }
            else
            {
                maxLabel <- c()
                maxLogLike <- -Inf
                for (j in 1:randomStart)
                {
                    label <- sample(1:K[i],ly,replace=T)
                    if (nu!=Inf)
                    {
                        obj<-.C("flowClust",
                            as.double(t(y)),
                            as.integer(ly),
                            as.integer(py),
                            as.integer(K[i]),
                            w=rep(0.0,K[i]),
                            mu=rep(0.0,K[i]*py),
                            precision=rep(0,K[i]*py*py),
                            lambda=as.double(lambda),
                            nu=as.double(nu),
                            z=rep(0.0,ly*K[i]),
                            u=rep(0.0,ly*K[i]),
                            as.integer(label),	
                            uncertainty=double(ly),
                            as.double(u.cutoff),
                            as.double(z.cutoff),
                            flagOutliers=integer(ly),
                            as.integer(B.init),
                            as.double(tol.init),
                            as.integer(trans),
                            logLike=as.double(0),
                            package="flowClust")
                    }
                    else
                    {
                        obj<-.C("flowClustGaussian",
                            as.double(t(y)),
                            as.integer(ly),
                            as.integer(py),
                            as.integer(K[i]),
                            w=rep(0.0,K[i]),
                            mu=rep(0.0,K[i]*py),
                            precision=rep(0,K[i]*py*py),
                            lambda=as.double(lambda),
                            z=rep(0.0,ly*K[i]),
                            u=rep(0.0,ly*K[i]),
                            as.integer(label),	
                            uncertainty=double(ly),
                            as.double(q.cutoff),
                            as.double(z.cutoff),
                            flagOutliers=integer(ly),
                            as.integer(B.init),
                            as.double(tol.init),
                            as.integer(trans),
                            logLike=as.double(0),
                            package="flowClust")
                    }
                    if (obj$logLike > maxLogLike)
                    {
                        maxLabel <- label
                        maxLogLike <- obj$logLike                    
                    }
                }
                label <- maxLabel
             }
        }
    
        if (nu!=Inf)
        {
            obj<-.C("flowClust",
                as.double(t(y)),
                as.integer(ly),
                as.integer(py),
                as.integer(K[i]),
                w=rep(0.0,K[i]),
                mu=rep(0.0,K[i]*py),
                precision=rep(0,K[i]*py*py),
                lambda=as.double(lambda),
                nu=as.double(nu),
                z=rep(0.0,ly*K[i]),
                u=rep(0.0,ly*K[i]),
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
        }
        else
        {
            obj<-.C("flowClustGaussian",
                as.double(t(y)),
                as.integer(ly),
                as.integer(py),
                as.integer(K[i]),
                w=rep(0.0,K[i]),
                mu=rep(0.0,K[i]*py),
                precision=rep(0,K[i]*py*py),
                lambda=as.double(lambda),
                z=rep(0.0,ly*K[i]),
                u=rep(0.0,ly*K[i]),
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
        sigma<-array(0,c(K[i],py,py))
        precision<-matrix(obj$precision,K[i],py*py,byrow=TRUE)
        for(k in 1:K[i]) sigma[k,,]<-matrix(precision[k,],py,py,byrow=TRUE)
    
        BIC <- 2*obj$logLike - log(ly)*(K[i]*(py+1)*py/2+K[i]*py+K[i]-1+as.integer(trans))
        z<-matrix(obj$z,ly,K[i],byrow=TRUE)
        ICL <- BIC + 2*sum(z*log(z),na.rm=TRUE)
    
        z <- u <- matrix(NA, length(include), K[i])
        z[include,] <- matrix(obj$z,ly,K[i],byrow=TRUE)
        u[include,] <- matrix(obj$u,ly,K[i],byrow=TRUE)
        tempLabel <- label
        label <- uncertainty <- flagOutliers <- rep(NA, length(include))
        label[include] <- tempLabel
        uncertainty[include] <- obj$uncertainty
        flagOutliers[include] <- as.logical(obj$flagOutliers)

        result[[i]] <- new("flowClust", expName=expName, varNames=varNames, K=K[i],
            w=obj$w, mu=matrix(obj$mu,K[i],py,byrow=TRUE), sigma=sigma, lambda=obj$lambda,
            nu=nu, z=z, u=u, label=label, uncertainty=uncertainty,
            ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min), 
            rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL)

    }
    
    if (length(K)==1)
    {
        result[[1]]
    }
    else
    {
        result <- new("flowClustList", result, criterion=criterion)    
        result@index <- which.max(criterion(result, criterion))
        result
    }    
    
}
