# 
#  mkPriors.R
#  flowClust with Bayes priors
#  
#  Created by Greg Finak on 2011-02-15.
#  Copyright 2011 Greg Finak. All rights reserved.
# 

#'Generate a prior specification based on a flowClust model
#
#' This function generates a prior specification based on a flowClust fit object
#' It can be passed to a second round of flowClust() with usePrior="yes" 
#' The prior could be estimated from a single sample, for example, and then
#' used to speed up the convergence for other samples. 
#'@param x a flowClust fit object
#'@param kappa is the fraction of equivalent observations by which to weight this prior relative to the flowClust model.
#'@param NT the number of total equivalent observation
#'@param addCluster not currently supported
#'@export
flowClust2Prior<-function(x,kappa,NT=10000,addCluster=NULL){
    Nt<-nrow(x@z)
    p<-ncol(x@mu)
    K<-x@K
    
    nu0<-Ng<-x@w*Nt
    if(all((nu0*kappa-p-1)>0)){
        Lambda0<-x@sigma;
        for(i in 1:K){
            Lambda0[i,,]<-Lambda0[i,,]*(kappa*nu0[i]-p-1)
        }
    }else{
        stop("Can't proceed. Prior nu0 is negative for cluster(s) ",paste(which((nu0-p-1)>0),collapse=","),"\n(p-1) = ",p-1,": Try increasing kappa")
    }
    
    Omega0<-array(0,c(K,p,p))
    for(i in 1:K){
        #Make precision of prior means spherical.
        Omega0[i,,]<-diag(1,p)
        if(p==1){
            dS<-x@sigma[i,,]
            dO<-Omega0[i,,]
        }else{
            dS<-det(x@sigma[i,,])
            dO<-det(Omega0[i,,])
        }
        k<-(dO/dS)^(1/p)
        #rescale Omega0
        Omega0[i,,]<-Omega0[i,,]*k
        #Multiply by Ng*kappa
        Omega0[i,,]<-solve(Omega0[i,,]*Ng[i]*kappa)
    }
    nu0<-nu0*kappa
    Mu0<-x@mu
    
    lambda<-x@lambda
    #w0 are dirichlet prior parameters
    #Make them vague
    w0<-x@w*Nt
    
    #if(!is.null(addCluster)){
    #        S<-cov(Mu0)
    #        Lam<-array(0,c(K+1,p,p))
    #        om<-array(0,c(K+1,p,p))
    #        Mu0<-rbind(Mu0,colMeans(Mu0))
    #        for(i in 1:K){
    #            om[i,,]<-Omega0[i,,]
    #            Lam[i,,]<-Lambda0[i,,]
    #        }
    #        om[K+1,,]<-diag(1,p)
    #        Lam[K+1,,]<-S
    #        if(p==1){
    #            dS<-Lam[K+1,,]
    #            dO<-om[K+1,,]
    #        }else{
    #            dS<-det(Lam[K+1,,])
    #            dO<-det(om[K+1,,])
    #        }
    #        k<-(dO/dS)^(1/p)
    #        om[K+1,,]<-om[K+1,,]*k
    #        om[K+1,,]<-solve(om[K+1,,])
    #        Omega0<-om
    #        Lambda0<-Lam
    #        nu0<-c(nu0,p+2)
    #        w0<-c(w0,1)
    #        K<-K+1
    #    }
    prior<-list(Mu0=Mu0,Lambda0=Lambda0,Omega0=Omega0,w0=w0,nu0=nu0,nu=x@nu,lambda=x@lambda,K=K)
    class(prior)<-"flowClustPrior"
    attr(prior,"lambda")<-x@lambda
    prior;
}

# =======================
# = Generic for mkPrior =
# =======================
setGeneric("mkPrior",function(gate,data,nu0,Omega0,...){
standardGeneric("mkPrior");	
})
# ===============================================================================
# = Generate a prior from a polygonGate and a flowFrame. Provide nu0 and Omega0 =
# ===============================================================================
setMethod("mkPrior",signature("polygonGate","flowFrame","numeric","matrix"),function(gate,data,nu0,Omega0){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);	
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-cov(data)*(nu0-length(dims)-1)
	Mu0<-colMeans(data)
	if(length(dims)!=ncol(Omega0)){
		stop("Invalid dimensions of \"Omega0\". Received ",dim(Omega0)[1]," x ",dim(Omega0)[2]," but expecting ",length(dims)," x ",length(dims));
	}
	prior=list(Mu0=Mu0,Lambda0=Lambda0,nu0=nu0,Omega0=Omega0)
	#message("Making prior from gate and data");
	prior;
})
#hyperparemeters not specified. Default nu0=500, Omega0=diag(1/1000,D)

# =========================================================
# = Nu0 and Omega0 specified. RectangleGate and flowFrame =
# =========================================================
setMethod("mkPrior",signature("rectangleGate","flowFrame","numeric","matrix"),function(gate,data,nu0,Omega0){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-cov(data)*(nu0-length(dims)-1)
	Mu0<-colMeans(data)
	if(length(dims)!=ncol(Omega0)){
		stop("Invalid dimensions of \"Omega0\". Received ",dim(Omega0)[1]," x ",dim(Omega0)[2]," but expecting ",length(dims)," x ",length(dims));
	}
	prior=list(Mu0=Mu0,Lambda0=Lambda0,nu0=nu0,Omega0=Omega0)
	prior;
})

# =========================================================================
# = hyperparemeters not specified. 
# = Returns an  incomplete prior specification. Shouldn't be called       =
# = Directly by the user.												  =
# =========================================================================
setMethod("mkPrior",signature("rectangleGate","flowFrame","missing","missing"),function(gate,data,nu0=NA,Omega0=NA){
	gc(reset=TRUE)
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-cov(data) #This is the covariance matrix, not the real Lambda0
	Mu0<-colMeans(data)
	n<-dim(data)[1];
	prior<-list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})
# ===================================================================================
# = Generate a prior from a polygonGate and a flowFrame with nu0 and Omega0 missing
# = Returns and incomplete prior specification. Should not be called by the user    =
# ===================================================================================
setMethod("mkPrior",signature("polygonGate","flowFrame","missing","missing"),function(gate,data,nu0=NA,Omega0=NA){
	##Extract the right dimensions
	dims<-unlist(lapply(gate@parameters,function(x)x@parameters),use.names=FALSE);
	dims<-dims[na.omit(match(parameters(data)@data$name,dims))]
	data<-exprs(data[,dims])
	Lambda0<-cov(data)
	Mu0<-colMeans(data)
	n<-dim(data)[1]
	prior=list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})


# =============================
# = S4 plot method for priors =
# =============================
#setMethod("plot",signature=c("flowClustPriorTree","GatingHierarchy"),function(x,y,node,dims=NULL,parent=TRUE,jitter=F,...){
#	if(is.numeric(node)){
#		if(node<=length(nodes(x))){
#			d<-nodeData(x,RBGL::tsort(x)[node],"gate")[[1]];
#			if(class(d)=="BooleanGateListList"){
#				message("node ", RBGL::tsort(x)[node], " is a Boolean Gate. There is no prior to plot.");
#			}else{
#				prior<-nodeData(x,RBGL::tsort(x)[node],"prior")[[1]];
#				if(parent){
#					data<-getData(y,getParent(y,getNodes(y,tsort=T)[node]))
#				}else{
#					data<-getData(y,node,tsort=TRUE);
#				}
#				plotPrior(data,prior,dims=dims);
#			}
#		}
#	}else{
#		stop("node must be numeric")
#	}
#})
# =============================
# = S4 plot method for priors =
# =============================
#setMethod("plot",signature=c("flowClustPriorTree","GatingSet"),function(x,y,node,dims=NULL,parent=TRUE,jitter=F,add=F,...){
#	if(is.numeric(node)){
#		if(node<=length(nodes(x))){
#			if(class(nodeData(x,RBGL::tsort(x)[node],"gate")[[1]])=="BooleanGateListList"){
#				message("node ", tsort(x)[node], " is a Boolean Gate. There is no prior to plot");
#				invisible(0);
#			}
#			prior<-nodeData(x,RBGL::tsort(x)[node],"prior")[[1]];
#			dims<-colnames(prior$Mu0)
#			data<-getData(y,node,tsort=TRUE)
#			if(!add){
#				if(parent){
#					pnode<-match(getParent(y[[1]],getNodes(y[[1]],tsort=T)[node]),getNodes(y[[1]],tsort=T))
#					if(jitter){
#							plot(jitter(fsApply(getData(y,pnode,tsort=T),function(x)exprs(x[,dims]))),pch='.',...);
#						}else{
#							plot(fsApply(getData(y,pnode,tsort=T),function(x)exprs(x[,dims])),pch='.',...);
#						}
#				}
#				else{
#					plot(fsApply(data,function(x)exprs(x[,dims])),pch='.',...)
#				}
#			}
#			l<-fsApply(data,function(x){
#					e<-exprs(x[,dims])
#					cm<-t(as.matrix(colMeans(e)))
#					el<-ellipse(cov((e)),centre=cm)
#					list(centers=cm,cvs=el)
#			})
#			lapply(l,function(x){points(x$centers,col="blue",pch="X",font=2);lines(x$cvs,col="gray",lwd=1)})
#			points(as.matrix(prior$Mu0),col="red",pch=20);
#			lines(ellipse(solve(prior$Omega0[1,,]),centre=prior$Mu0),col="red",lwd=1,lty=2);
#			lines(ellipse(prior$Lambda0[1,,]/(prior$nu0-2-1),centre=prior$Mu0),col="red",lwd=1,lty=1)
#		}
#	}else{
#		stop("node must be numeric")
#	}
#})

# =========================================================================
# = Plot a prior given some data (a flowFrame) and a prior specification. =
# = The prior specification is for a single cluster. 					  =
# =========================================================================

plotPrior<-function(data,prior,dims=NULL,...){
	if(!"smooth"%in%names(list(...))){
		sm=FALSE
	}
	if(class(prior)=="flowClustPriorList"){
		prior<-prior[[1]];
	}
	if(!is.null(dims)){
		if(all(dims%in%colnames(prior$Mu0))){
			dims<-dims[na.omit(match(colnames(prior$Mu0),dims))]
			dim.inds<-match(dims,colnames(prior$Mu0))
		}else{
			stop("Can't find ",dims[which(!(dims%in%colnames(prior$Mu0)))]," in the prior.");
		}
	}else{
		dims<-colnames(prior$Mu0)
		dim.inds<-match(dims,colnames(prior$Mu0))
	}
	k<-nrow(prior$Mu0);
	nd<-ncol(prior$Mu0)
	if(nd>1){
		if (exists("sm")){
			flowViz:::fplot(data[,dims],smooth=sm,...);
		}else{
			flowViz:::fplot(data[,dims],...);
		}			
		for(i in 1:k){
			points(t(as.matrix(prior$Mu0[i,dim.inds])),pch=20,col="red")
			lines(ellipse(solve(prior$Omega0[i,dim.inds,dim.inds]),centre=prior$Mu0[i,dim.inds]),col="green",lwd=2,lty=2)
			lines(ellipse(prior$Lambda0[i,dim.inds,dim.inds]/(prior$nu0[i]-nd-1),centre=prior$Mu0[i,dim.inds]),col="red",lwd=2,lty=2)
		}
	}else{
		for(i in 1:k){
			if (exists("sm")){
				flowViz:::fplot(data[,dims],smooth=sm,breaks=256,...);
			}else{
				flowViz:::fplot(data[,dims],breaks=256,...);
			}
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds])),col="red")
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.975)*sqrt(1/prior$Omega0[i,dim.inds,dim.inds]),col="red",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.025)*sqrt(1/prior$Omega0[i,dim.inds,dim.inds]),col="red",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.975)*sqrt((prior$Lambda0[i,dim.inds,dim.inds])/(prior$nu0-2)),col="green",lty=2);
			abline(v=t(as.matrix(prior$Mu0[i,dim.inds]))+qnorm(0.025)*sqrt((prior$Lambda0[i,dim.inds,dim.inds])/(prior$nu0-2)),col="green",lty=2);
		}
	}
	invisible(0);
}
# ===========================================================================================
# = Pad the prior to K clusters. K must be greater than the number of clusters in the prior =
# ===========================================================================================
.padPriorToK<-function(prior,k,env){
	#env is an environment with the data to be used to pad the prior.
	flag<-0;
	if(inherits(prior,"flowClustPriorList")){
		prior<-prior[[1]];
		flag<-1;
	}else if(inherits(prior,"flowClustPrior")){
		prior<-prior;
	}else{
		stop("prior must be of class \"flowClustPrior\" or \"flowClustPriorList\"");
	}
	nd<-ncol(prior$Mu0);
	nc<-nrow(prior$Mu0);
	if(nc>=k){
		stop("Prior cannot be padded to fewer than or equal to ", k, " clusters. It already has ",nc," clusters.");
	}
	Mu0<-matrix(NA,k,nd);
	Omega0<-array(NA,c(k,nd,nd));
	Lambda0<-array(NA,c(k,nd,nd));
	nu0<-numeric(k);
	colnames(Mu0)<-colnames(prior$Mu0);
	for(i in 1:nc){
		Mu0[i,]<-prior$Mu0[i,];
		Omega0[i,,]<-prior$Omega0[i,,]
		Lambda0[i,,]<-prior$Lambda0[i,,]
		nu0[i]<-prior$nu0[i]
	}
	for(i in (nc+1):k){
		Mu0[i,]<-rep(0,nd)##Set to vague values
		Omega0[i,,]<-diag(0,nd)## Set to vague values
		Lambda0[i,,]<-diag(0,nd)## Set to vague values
		nu0[i]<-nd+1; #vague
	}
	pprior<-list(Mu0=Mu0,Omega0=Omega0,Lambda0=Lambda0,nu0=nu0);
	class(pprior)<-"flowClustPrior";
	if(flag){
		pprior<-list(pprior);
		class(pprior)<-"flowClustPriorList";
	}
	return(pprior);
}

#setMethod("getParent",signature=c("flowClustTree","character"),function(obj,y){
#	flowClust:::.graph_parent(obj,y)
#})
# ========================================================
# = Make a tree of prior distributions from a Gating Set =
# ========================================================
# TODO Make the returned object a specific class (priorTree)
#mkPriorTree<-function(gh,model.cov="full",...){
#	tree<-gh[[1]]@tree;
#	q<-RBGL::tsort(tree)[1]
#	allnodes<-RBGL::tsort(tree)
#	g<-new("flowClustPriorTree")
#	while(length(q)!=0){
#		#start at the root
#		current<-q[1];
#		q<-q[-(1)];
#		curind<-na.omit(match(current,allnodes));
#		g<-graph::addNode(current,g)
#		parent<-.graph_parent(tree,current)
#		
#		#enqueue all children
#		q<-c(q,unlist(adj(tree,current)));
#		if(length(parent)!=0){
#			g<-graph::addEdge(parent,current,g);
#		}
#		#Does the current node have a gate?
#		if(all(is.na(flowWorkspace::getGate(gh,curind,tsort=TRUE)))){
#			#No?
##			 #message("node ",current," has no gate")			
#			next;
#		}else{
#			#Yes?
#			#message("node ",current, " has a gate");
#			gate<-flowWorkspace::getGate(gh,curind,tsort=TRUE);
#			if(all(unlist(lapply(gate,function(x)class(x)=="BooleanGate")))){
#				#message("Omitting ",gate[[1]]$ref);
#				class(gate)<-"BooleanGateList"
#				gate<-list(gate);
#				class(gate)<-"BooleanGateListList"
#				# Priors are now computed in 2D on the parent projection.. clearly this is not ideal. Should use the 
#				# boolean gate defining dimensions.
#				#Perhaps it would be better to take a boolean combination of the regular gates, rather than use a prior.
#				#dims<-getDimensions(gh[[1]], tsort(gh[[1]]@tree)[curind])
#				#dim.ind<-match(selectMethod("colnames",signature("flowFrame"))((getData(gh[[1]], curind, tsort = TRUE))),dims)
#				#dims<-dims[na.omit(dim.ind)]
#				#prior<-mkPrior(data=getData(gh,curind,tsort=TRUE)[,dims],nu0=NA,model.cov="full");
#			}else{
#				whp<-which(flowWorkspace::getNodes(gh[[1]],tsort=TRUE)%in%.graph_parent(gh[[1]]@tree,flowWorkspace::getNodes(gh[[1]],tsort=TRUE)[curind]))
#				prior<-mkPrior(gate,gh,nu0=NA,model.cov=model.cov,indexp=whp,indexc=curind)	
#				gc(reset=T)				
#				nodeData(g,current,"prior")<-prior;
#			}
#			gate<-list(gate)
#			class(gate)<-"GateList"
#			nodeData(g,current,"gate")<-gate;
#		}
#		gc(reset=TRUE)
#	
#	}
#	cat("\nDone!\n");
#	g
#}


# ============================================================
# = Get the parent node of a given node in a graphNEL object =
# ============================================================
.graph_parent<-function(g,n){
	return(setdiff(unlist(graph::adj(ugraph(g),n)),unlist(graph::adj(g,n))))
}

# ==================================================================================================
# = Estimate the hyperparameters of the prior given the means and covariances of multiple samples. =
# ==================================================================================================
.estimateHyperParameters<-function(priors,model.means,model.cov,nu0){
	#Empirical Bayes Estimation of hyperparameters (except nu0)
	d<-dim(priors[[1]]$Lambda0)[1]
	MuG0<-NA;LambdaG0<-NA;OmegaG0<-NA;nuG0<-NA;
	MuG0<-colMeans(do.call(rbind,lapply(priors,function(x)x$Mu0)))
	if(model.means=="full"){
		#What if n<p
		cm<-do.call(rbind,lapply(priors,function(x)x$Mu0))
		if(nrow(cm)==1){
			
			# TODO Test this more thoroughly, i.e. case where n=1, and using covariance of the means at 1% of the covariance of the data.
			OmegaG0<-solve(diag(diag(priors[[1]]$Lambda*0.01)))
		}else if(nrow(cm)<ncol(cm)){
			OmegaG0<-solve(cov.shrink(cm,verbose=FALSE));
		}else{
			OmegaG0<-try(solve(cov(cm)),silent=TRUE);
			if(inherits(OmegaG0,"try-error")){
				OmegaG0<-solve(cov.shrink(cm,verbose=FALSE))
			}
		}
		#Sometimes the estimate is unstable and negative. In that case use a diagonal covariance parameterization
		if(OmegaG0[1,1]<0){
			OmegaG0<-solve(diag(apply(cm,1,var)))
		}
		
	
	}else if(model.means=="DU"){
		# TODO code to handle nrow = 1
		d<-dim(priors[[1]]$Lambda0)[1]
		OmegaG0<-solve(diag(diag(cov(do.call(rbind,lapply(priors,function(x)x$Mu0)))),d))
#		if(OmegaG0[1,1]<0){
#		}
	}else if(model.means=="DE"){
		# TODO code to handle nrow=1
		d<-dim(priors[[1]]$Lambda0)[1]
		OmegaG0<-solve(diag((det(cov(do.call(rbind,lapply(priors,function(x)x$Mu0)))))^(1/d),d))
	}
	nuG0<-nu0	
	if(model.cov=="full"){
		#This estimate of nu0 is not correct for this model, working on it.
		if(is.na(nu0)){
			X<-lapply(priors,function(x)x$Lambda0)
			d<-dim(priors[[1]]$Lambda0)[1]		
			nu0<-optimize(f=.LIW,interval=c(d+2,20000),X=X)$minimum
			nuG0<-nu0
		}
		LambdaG0<-Reduce("+",lapply(priors,function(x)x$Lambda0))/length(priors)
		LambdaG0<-LambdaG0*(nu0) #ML estimate of LambdaG0
	}else if(model.cov=="DE"){
		#Estimate nu0
		if(is.na(nu0)){
			d<-dim(priors[[1]]$Lambda0)[1]
	nu0<-(solve(diag(diag(Reduce("+",lapply(priors,function(x)(x$Lambda0)))),d))*d)/(sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0)))))*(d-1))+diag(1/(d-1),d)
			nu0<-det(nu0)^(-1/d)
			nuG0<-nu0
		}
		#Or use the given value	
			d<-dim(priors[[1]]$Lambda0)[1]		
			LambdaG0<-diag((length(priors)*dim(priors[[1]]$Lambda0)[1]*nu0)/sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0))))),d)
	}else if(model.cov=="DU"){
		#This estimate of nu0 is not correct for this model.. working on it.
			if(is.na(nu0)){
				d<-dim(priors[[1]]$Lambda0)[1]	
				X<-lapply(priors,function(x)x$Lambda0)	#nu0<-(solve(diag(diag(Reduce("+",lapply(priors,function(x)(x$Lambda0)))),d))*d)/(sum(diag(Reduce("+",lapply(priors,function(x)solve(x$Lambda0)))))*(d-1))+1/(d-1)
				# nu0<-det(nu0)^(-1/d)
				nu0<-optimize(f=.LIW,interval=c(d+2,20000),X=X)$minimum
				nuG0<-nu0
			}
			d<-dim(priors[[1]]$Lambda0)[1]	
		LambdaG0<-diag(diag(Reduce("+",lapply(priors,function(x)x$Lambda0))/length(priors)),dim(priors[[1]]$Lambda0)[1])*(nuG0-d-1)
	}
	n<-list(unlist(lapply(priors,function(x)x$n)))# number of events in each sample. To be used to compute wi later.
	#include the number of events for dirichlet
	prior=list(Mu0=MuG0,Omega0=OmegaG0,Lambda0=LambdaG0,nu0=nuG0,n=n);
	# Mu0 must be a matrix for flowClust
	prior$Mu0<-t(as.matrix(prior$Mu0)) 
	#Omega0 and Lambda0 must be K,py,py arrays
	prior$Lambda0<-array(prior$Lambda0,c(1,ncol(prior$Lambda0),ncol(prior$Lambda0)))
	prior$Omega0<-array(prior$Omega0,c(1,ncol(prior$Omega0),ncol(prior$Omega0)))
	class(prior)<-"flowClustPrior";
	prior<-list(prior);
	class(prior)<-"flowClustPriorList";
	return(prior);
}

# ==========================================
# = Mulitple gates (list) and a gatingset  =
# ==========================================
#setMethod("mkPrior",signature("list","GatingSet",nu0="ANY","missing"),function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full",indexp=NULL,indexc=NULL){
#	if(is.null(indexp)|is.null(indexc)){
#		stop("Must specify for which node you want to construct a prior via indexp or indexc=\"numeric\" argument");
#	}
#	priors<-list();
#	method=match.arg(model.cov,c("full","DE","DU"));
#	model=match(model.means,c("full","DU","DE"))
#	repflag<-FALSE;
#	estflag<-FALSE;
#	if(!all(unlist(lapply(gate,function(x)class(x)=="polygonGate"|class(x)=="rectangleGate"),use.names=FALSE))){
#		#browser()
#		stop("mkPrior: All elements of \"gate\" must be class polygonGate or rectangleGate. Are the gating hierarchies in the gating set identical?");
#	}
#	if(!(is(nu0,"numeric")|is.na(nu0))){
#		stop("nu0 must be a numeric vector of length >= 1, or NA for estimation");
#	}
#	if(length(nu0)>1){
#		repflag<-TRUE;
#	}
#	cat(".")
#	#Deal only with gates where there are more than 5 samples
#	for(i in 1:length(data)){
#		if(nrow(getData(data[[i]],indexc,tsort=T))>5){
#			priors[[i]]<-mkPrior(gate[[i]],getData(data[[i]],indexc,tsort=TRUE))
#		}else{
#			#Ignore the gate if there's not enough sample points
#			next;
#		}
#	}
#	priors<-priors[unlist(lapply(priors,function(x)!is.null(x)))]
#	if(length(priors)==0){
#		#No gate has enough data to construct a prior.
#		return(NA);
#	}
#	prior<-.estimateHyperParameters(priors,model.means=model.means,model.cov=model.cov,nu0=nu0);
#	return(prior);
#})
# ======================================================================================
# = We should also be able to generate a prior from a single gate and multiple samples =
# ======================================================================================
 setMethod("mkPrior",signature("list","flowSet",nu0="missing","missing"),function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
 mkPrior(gate,data,nu0=NA)
 })

# ==============================
# = Prior from a flowSet alone =
# ==============================
setMethod("mkPrior",signature("missing","flowSet",nu0="ANY","missing"),function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
	priors<-list();
	method=match.arg(model.cov,c("full","DE","DU"));
	model=match(model.means,c("full","DU","DE"))
	repflag<-FALSE;
	estflag<-FALSE;
	if(!(is(nu0,"numeric")|is.na(nu0))){
		stop("nu0 must be a numeric vector of length >= 1, or NA for estimation");
	}
	if(length(nu0)>1){
		repflag<-TRUE;
	}
	cat(".")
	for(i in 1:length(data)){
		if(nrow(data[[i]])>5){
			priors[[i]]<-mkPrior(data=data[[i]])
		}else{
			#Ignore the gate if there's not enough sample points
			next;
		}
	}
	priors<-priors[unlist(lapply(priors,function(x)!is.null(x)))]
	if(length(priors)==0){
		#No gate has enough data to construct a prior.
		return(NA);
	}

	prior<-.estimateHyperParameters(priors,model.means=model.means,model.cov=model.cov,nu0=nu0);
	return(prior);
	
})
# ================================================================================
# = Construct a prior from a flowFrame alone. Not meant to be called by the user =
# ================================================================================
setMethod("mkPrior",signature("missing","flowFrame",nu0="missing","missing"),function(gate,data,nu0,Omega0){
	gc(reset=TRUE)
	##Use all the dimensions, since they're not specified.
	data<-exprs(data)
	if(ncol(data)>=nrow(data)){
		Lambda0<-cov.shrink(data,verbose=FALSE)
		class(Lambda0)<-"matrix"
	}else{
		#Lambda0<-cov.shrink(data,verbose=FALSE) #This is the covariance matrix, not the real Lambda0
		#class(Lambda0)<-"matrix"
		Lambda0<-cov(data) #This is the covariance matrix, not the real Lambda0
		
	}
	Mu0<-colMeans(data)
	n<-dim(data)[1];
	prior<-list(Mu0=Mu0,Lambda0=Lambda0,n=n)
	gc(reset=T)
	prior;
})
# ====================================================
# = Multiple gates (same gates) and multiple samples =
# = Calls mkPrior with the "missing" hyperparamter   =
# = Signatures										 =
# ====================================================
setMethod("mkPrior",signature("list","flowSet",nu0="ANY","missing"),function(gate,data,nu0=NA,Omega0,model.cov="full",model.means="full"){
	priors<-list();
	method=match.arg(model.cov,c("full","DE","DU"));
	model=match(model.means,c("full","DU","DE"))
	repflag<-FALSE;
	estflag<-FALSE;
	if(!all(unlist(lapply(gate,function(x)class(x)=="polygonGate"|class(x)=="rectangleGate"),use.names=FALSE))){
		stop("All elements of \"gate\" must be class polygonGate or rectangleGate");
	}
	if(!(is(nu0,"numeric")|is.na(nu0))){
		stop("nu0 must be a numeric vector of length >= 1, or NA for estimation");
	}
	if(length(nu0)>1){
		repflag<-TRUE;
	}
	cat(".")
	#Deal only with gates where there are more than 5 samples
	for(i in 1:length(data)){
		sub<-Subset(data[[i]],filter(data[[i]],gate[[i]]))
		if(nrow(sub)>5){
			priors[[i]]<-mkPrior(gate[[i]],sub)
		}else{
			#Ignore the gate if there's not enough sample points
			next;
		}
	}
	priors<-priors[unlist(lapply(priors,function(x)!is.null(x)))]
	if(length(priors)==0){
		#No gate has enough data to construct a prior.
		return(NA);
	}

	prior<-.estimateHyperParameters(priors,model.means=model.means,model.cov=model.cov,nu0=nu0);
	return(prior);
})

# ====================================
# = Log-likelihood for the dirichlet =
# ====================================
.LD<-function(alpha,x){
	
	-sum(log(MCMCpack::ddirichlet(x,alpha)))
}

# =======================================================
# = Log-Likelihood for the Inverse Wishart Distribution =
# =======================================================
.LIW<-function (nu, X) 
{
    n <- length(X)
    d <- nrow(X[[1]])
    Sinv <- lapply(X, function(x) t(solve(x)))
    sumSinv <- Reduce("+", lapply(Sinv, function(x) t(x)))
    S = solve(sumSinv/(nu * n))
    A <- n * d * (d - 1)/4 * log(pi) - n * sum((sapply(1:d, function(j) lgamma(nu/2 + 
        (1 - j)/2))))
    B <- n * nu/2 * log(det(S))
    C <- -Reduce("+", lapply(X, function(x) log(det(x)))) * (nu - 
        d - 1)/2
    D <- -Reduce("+", lapply(Sinv, function(x) sum(diag(S %*% 
        x))/2))
    E <- -n * nu * d/2 * log(2)
    return(-sum(c(A, B, C, D, E)))
}
# ======================
# = Combine two priors =
# ======================
.mergePriors<-function(priors){
	l<-length(priors)
	nd<-dim(priors[[1]]$Mu0)[2]
	Omega0<-array(NA,c(l,nd,nd));
	Lambda0<-array(NA,c(l,nd,nd));
	nu0<-numeric(l);
	Mu0<-matrix(NA,l,nd);
	w0<-numeric(l);
	for(i in 1:l){
		Omega0[i,,]<-priors[[i]]$Omega0[1,,];
		Lambda0[i,,]<-priors[[i]]$Lambda0[1,,];
		Mu0[i,]<-priors[[i]]$Mu0[1,];
		nu0[i]<-priors[[i]]$nu0[1];	
	}
	#estimate w0 from dirichlet likelihood. Numerical optimization. Initialized at unity.
	x<-do.call(cbind,unlist(lapply(priors,function(x)x$n),recursive=F))
	for(i in 1:nrow(x)){
		x[i,]<-x[i,]/sum(x[i,])
	}
	wstar<-optim(par=rep(1,length(w0)),x=x,fn=.LD,method="L-BFGS-B",lower=rep(1,length(w0)))
	if(wstar$convergence!=0){
		w0<-rep(1,length(w0));
	}
	else{
		w0<-wstar$par;
	}
	
	colnames(Mu0)<-colnames(priors[[1]]$Mu0)
	mprior<-list(Mu0=Mu0,Omega0=Omega0,Lambda0=Lambda0,nu0=nu0,w0=w0);
	class(mprior)<-"flowClustPrior";
	return(mprior)
}
# ============================================================================================================
# = Determine if the priors for the children of a given node should be combined into a single specification. =
# ============================================================================================================
.combinePriorsForChildNodes<-function(x,n){
	children<-unlist(adj(x,n),use.names=F);
	priors<-nodeData(x,children,"prior");
	
	#Which children have a  prior specification?
	havePriors<-unlist(lapply(priors,function(x)inherits(x,"flowClustPrior")))
	#Compare only those with priors
	priors<-priors[havePriors];
	if(length(priors)>0){
	dims<-lapply(priors,function(x)colnames(x$Mu0))
	groups<-matrix(0,nrow=length(dims),ncol=length(dims))
	colnames(groups)<-names(dims);
	rownames(groups)<-names(dims);
	unassigned<-names(dims);
	assigned<-list();
	while(length(unassigned)!=0){
		i<-unassigned[1];
		k<-length(assigned)+1
		assigned[[k]]<-i;
		unassigned<-setdiff(unassigned,i)
		for(j in setdiff(unassigned,i)){
			if(identical(dims[[i]],dims[[j]])){
				assigned[[k]]<-c(assigned[[k]],j)
				unassigned<-setdiff(unassigned,j)
			}
		}		
	}
	mergedPriors<-list()
	#assigned is a list containing the names of identical nodes.
	## TODO fix this it doesn't appear to be correct ordering
	for(i in 1:length(assigned)){
		mergedPriors[[i]]<-.mergePriors(priors[assigned[[i]]])
		attr(mergedPriors[[i]],"popnames")<-assigned[[i]];
	}
	class(mergedPriors)<-"flowClustPriorList"
	}else{
		mergedPriors<-NA
	}
	mergedPriors;
}

#setGeneric("fitPriorTree2Sample",function(priortree,flowframe,...){
#	standardGeneric("fitPriorTree2Sample")
#})
#setMethod("fitPriorTree2Sample",signature=c("flowClustPriorTree","flowFrame"),function(priortree,flowframe,trans=0,nu.est=1,w0=c(5,1,0.5,0.01),nu0=NULL,rare.thresh=c(0.02,0.15,0.35),pi0=NULL,...){
#	.fitPriorTree2Sample(priortree=priortree,flowframe=flowframe,trans=trans,nu.est=nu.est,w0=w0,nu0=nu0,rare.thresh=rare.thresh,pi0=pi0,...)
#})

#setMethod("fitPriorTree2Sample",signature=c("flowClustPriorTree","GatingHierarchy"),function(priortree,flowframe,trans=0,nu.est=1,w0=c(5,1,0.5,0.01),nu0=NULL,rare.thresh=c(0.02,0.15,0.35),pi0=NULL,...){
#	.fitPriorTree2Sample(priortree=priortree,flowframe=flowframe,trans=trans,nu.est=nu.est,w0=w0,nu0=nu0,rare.thresh=rare.thresh,pi0=pi0,...)
#})

#.fitPriorTree2Sample<-function(priortree,flowframe,trans=0,nu.est=1,w0=c(5,1,0.5,0.01),nu0=NULL,rare.thresh=c(0.02,0.15,0.35),pi0=NULL,...){
	# if(is.null(names(add.nu0))&length(add.nu0)==1){
	# 	add.nu0<-rep(add.nu0,length(nodes(priortree)))
	# 	names(add.nu0)<-nodes(priortree);
	# }
	# if(is.null(names(scale.Omega0))&length(scale.Omega0)==1){
	# 	scale.Omega0<-rep(scale.Omega0,length(nodes(priortree)))
	# 	names(scale.Omega0)<-nodes(priortree);
	# }
	# if(is.null(names(scale.Lambda0))&length(scale.Lambda0)==1){
	# 	scale.Lambda0<-rep(scale.Lambda0,length(nodes(priortree)))
	# 	names(scale.Lambda0)<-nodes(priortree)
	# }
	#test if it's a gatinghierarchy.. so we can fit to gated samples.
#	if(class(flowframe)=="GatingHierarchy"){
#		flowframe<-getData(flowframe);
#	}
#	fct<-new("flowClustTree",priortree);
	#A queue containing the root;
#	q<-RBGL::tsort(priortree)[1];
	#Add the root to the tree
#	fct<-graph::addNode(q,fct);
	#Add the data to the root
#	nodeData(fct,q[1],"population")<-flowframe;
#	while(length(q)!=0){
#		X<-q[1];
		#parent<-.graph_parent(priortree,X)
#		q<-q[-1L];
		#Get child nodes and use them to name the populations.
#		children<-unlist(adj(priortree,X),use.names=FALSE)
		#The current node's data.
#		flowframe<-nodeData(fct,X,"population")[[1]]
#		if(length(children)!=0){
			#enqueue the children
#			q<-c(q,children);
			
			#get the priors for the children
			# this is okay, since combinePriors will omit nodes where a prior is not defined.
#			priors<-flowClust:::.combinePriorsForChildNodes(priortree,X)
			#Fit flowClust to all the children (presumably the current node was fitted at the last iteration)
			#check if the flowframe is not empty and if there are child nodes with priors (ie not just boolean gates)
#			if(nrow(flowframe)!=0&!all(is.na(priors))){
				#nextpnames<-unlist(lapply(priors,function(x)attr(x,"popnames")))
				#res<-flowClustBayes(x=flowframe,prior=priors,trans=trans,nu.est=nu.est,add.nu0=add.nu0[nextpnames],scale.Lambda0=scale.Lambda0[nextpnames],scale.Omega0=scale.Omega0[nextpnames],...);
#				res<-flowClustBayes(x=flowframe,prior=priors,trans=trans,nu.est=nu.est,w0=w0,nu0=nu0,rare.thresh=rare.thresh,pi0=pi0,...);
				#lapply(res,function(x)cat("gating ",attr(x@prior,"popnames"),"\n"))
				# TODO This reordering of names should occur within the flowClust code.
#				res<-lapply(res,function(r){attr(r@prior,"popnames")<-attr(r@prior,"popnames")[r@prior$order];r})
				#res is of length length(priors);
#				pops<-lapply(res,function(r){
					# How should we represent an empty population?
					# just use the code below, omit the split method which throws errors when a population is missing.
					#pops<-try(selectMethod("split",c("flowFrame","flowClust"))(flowframe,r))
					#if(inherits(pops,"try-error")){
#						pops<-vector("list",r@K)
#						for(kk in 1:r@K){
							#browser();
#							pops[[kk]]<-flowframe[which(Map(r)==kk),]
#						}
					#}
#						names(pops)<-attr(r@prior,"popnames");
#						pops
#				})
#			}else{
#				pops<-vector("list",length(children))
#				pops<-sapply(pops,function(x)flowframe)
#				names(pops)<-children
#				pops<-list(pops)
#				res<-list();
#				
#			}
			#add children
#			fct<-graph::addNode(children,fct)
			#Add edges to the current node
#			fct<-graph::addEdge(X,children,fct)
#		
			#fill in the model for each child node.
#			if(length(res)!=0){
#				for(i in 1:length(res)){
#					e<-new.env(parent=emptyenv())
#					for(j in attr(res[[i]]@prior,"popnames")){
#						nodeData(fct,j,"model")<-e
#					}
#					assign("model",res[[i]],e);
#				}
#			}
			#fill in the data for each child node
#			for(i in 1:length(pops)){
#				for(j in 1:length(pops[[i]])){
#					nodeData(fct,names(pops[[i]])[[j]],"population")<-pops[[i]][[j]]
#				}
#			}
			#Fill in names
#		  	for(i in 1:length(pops)){
#				for(j in 1:length(pops[[i]])){
#					nodeData(fct,names(pops[[i]])[[j]],"name")<-names(pops[[i]])[[j]]
#				}
#		  	}
			#Get the gate for each node and add to graph.
#			for (x in children){
#				g<-nodeData(priortree,x,"gate")[[1]];
#				if(class(g)=="BooleanGateListList"){
					#Add the definition of the first gate, since the node names will match this tree.
					#Other gates are identical, but node names are different. Don't need them.
#					nodeData(fct,x,"gate")<-g[[1]][1];
#				}else{
#					nodeData(fct,x,"gate")<-g[1];
#				}
#			}
#		}
		# define populations for boolean gates below the current node
		# Add gates from the prior tree  to the flowclust tree.
		#These are the boolean nodes. They would have been skipped by flowClust above.
#		bnodes<-children[unlist(lapply(nodeData(priortree)[children],function(x)class(x$gate)=="BooleanGateListList"),use.names=TRUE)]
		#for each boolean node, pull the gate definition and construct the data
		#This code constructs  
		#a string to be evaluated and returns the indices of the cells in the boolean gate
#		if(length(bnodes)!=0){
#		for(i in 1:length(bnodes)){
			#cat("gating ",bnodes[i],"\n");
			# TODO modify so that we can support boolean nodes that are combinations of other boolean nodes. 
			# TODO modify for combinations like A/B/C & A/B/X/D, so that we have indices of the same length based on node B.
			#gate definition
#			bg<-nodeData(fct,bnodes,"gate")[[i]];
			#common ancestor
			#TODO check if all ref are in fct.
			
#			sps<-sapply(bg$ref,function(x){
#				RBGL::sp.between(fct,nodes(fct)[1],x)
#			})
#			ca.ind<-min(unlist(lapply(sps,function(x)x$length+1)))
#			ca<-sps[[1]]$path_detail[min(unlist(lapply(sps,function(x)x$length+1)))]
#			amodel<-getModel(fct,ca)
#			indices<-list();
#			for(j in 1:length(sps)){
#				indices[[j]]<-rep(TRUE,dim(amodel@z)[1])
				#indices[[j]]<-indices[[j]][which(indices[[j]]==T)]
#				for(k in (ca.ind):(sps[[j]]$length+1)){ 
#					this<-sps[[j]]$path_detail[k]
#					
#					thismodel<-getModel(fct,this);
#					if(!is.null(thismodel)){
#						ti<-match(this,attr(getModel(fct,this)@prior,"popnames"));
#						indices[[j]][which(indices[[j]]==TRUE)]<-indices[[j]][which(indices[[j]]==TRUE)]&Map(getModel(fct,this))==ti
#					}else{
#						#model at the current node is null, therefore the previous flowframe is probably empty.
#						#None of the events are part of this flowframe
#						indices[[j]][!is.na(indices[[j]])]<-FALSE
#					}
#				}	
#			}
#			query<-vector("list",length(indices));
#			query<-lapply(query,function(x)character())
#			for(j in 1:length(bg$ref)){
#				query[[j]]<-paste("indices[[",j,"]]",sep="")
#				query[[j]]<-paste(bg$v[j],query[[j]],sep="")
#			}
#			for(j in 1:length(bg$v2)){
#				query[[j+1]]<-(paste(unlist(query)[j:(j+1)],collapse=bg$v2[j]))
#			}
#			query<-query[[length(query)]];
#			#These are the indices of events included in the boolean gate relative to the parent of the component gates.
#			indices<-eval(parse(text=query))
#			indices[is.na(indices)]<-FALSE
#			d<-getData(fct,ca)[indices,]
#			nodeData(fct,bnodes[i],"population")<-d;
#		}
#		}
#		
#	}
#	return(fct);
#}
# ======================================================================
# = S4 method for plotting populations below a node in a flowClustTree =
# ======================================================================
#TODO correct axis for raw scale, rather than transformed scale.
#setMethod("plot",c("flowClustTree","character"),function(x,y,level=0.9,pch='.',...){
#	p<-flowClust:::.graph_parent(x,y);
#	if(class(nodeData(x,y,"gate")[[1]])=="BooleanGate"&length(p)!=0){
#		#plot the points in the gate using the projection of the parent
#		 selectMethod("plot",c("flowFrame","missing"))(getData(x,y,parent=F)[,getModel(x,p)@varNames],col="red",pch='x',smooth=FALSE,...)
#	}else{
#		if(length(p)==0){
#			selectMethod("plot",c("flowFrame","missing"))(getData(x,y,parent=FALSE),...)
#		}else{
#			#Need the range
#			data<-getData(x,y);
#			model<-getModel(x,y);
#			if(is.null(model)){
#				# #plot based on the prior.
#				# prior<-nodeData(x,y,"prior")[[1]]
#				# varnames<-colnames(prior$Mu0);
#				# d<-ncol(prior$Mu0);
#				# K<-nrow(prior$Mu0);
#				# range<-range(data[,varnames])
#				# plot(prior$Mu0,xlim=range[,varnames[1]],ylim=range[,varnames[2]],pch="x");
#				# sapply(1:K,function(i)lines(ellipse(prior$Omega0[,,i]),lwd=2,lty=2));
#				# sapply(1:K,function(i)lines(ellipse(prior$Lambda0[,,i]/prior$nu0[i]),lwd=2,lty=1));
#				frame();
#			}else{
#				range<-range(data[,model@varNames])
#				#parameter names
#				xlim<-range[,model@varNames[1]];
#				xlab<-pData(parameters(data))$desc
#				names(xlab)<-pData(parameters(data))$name
#				xlab<-xlab[model@varNames[1]]			
#				#dealing with missing values in parameters names
#				if(is.na(xlab)){
#					xlab<-model@varNames[1]
#				}else if(xlab==" "|xlab==""){
#					xlab<-model@varNames[1]
#				}
#				if(length(model@varNames)==2){
#					ylim<-range[,model@varNames[2]]
#					ylab<-pData(parameters(data))$desc
#					names(ylab)<-pData(parameters(data))$name
#					ylab<-ylab[model@varNames[2]]
#					if(is.na(ylab)){
#						ylab<-model@varNames[2]
#					}else if(ylab==" "|ylab==""){
#						ylab<-model@varNames[2]
#					}
#					selectMethod("plot",c("flowClust","missing"))(model,data=data,level=level,pch=pch,xlim=xlim,ylim=ylim,npoints=10000,xlab=xlab,ylab=ylab,...)
#				}else{
#					selectMethod("plot",c("flowClust","missing"))(model,data=data,level=level,pch=pch,xlim=xlim,npoints=10000,xlab=xlab,...)
#				}
#			}
#		}
#	}
#})

#setModel<-function(x,y,z){
#	assign("model",z,envir=nodeData(x,y,"model")[[1]])
#}
# TODO rewrite so that we can set the level argument for the model (quantile)

#setMethod("getPopStats",signature=c("flowClustTree"),function(x,level=0.9){
#	d<-t(sapply(bfs(x),function(y){
#		m<-getModel(x,y);
#	#	ruleOutliers(m)<-list(level=level);
#	#	setModel(x,y,m);
#		if(inherits(m,"flowClust")){
#			wh<-which(attr(m@prior,"popnames")%in%y);
#num<-length(which(Map(m)==wh))	
#den<-length(Map(m))	
#c(name=as.vector(y),model.proportion=m@w[wh],proportion=sapply(num/den,function(q)ifelse(is.infinite(q)|is.nan(q),0,q)),parent=den,count=num)
#		}else{
#		num<-nrow(getData(x,y,parent=FALSE))
#		den<-nrow(getData(x,y,parent=TRUE))	
#		c(name=as.vector(y),model.proportion=sapply(num/den,function(q)ifelse(is.infinite(q)|is.nan(q),0,q)),proportion=sapply(num/den,function(q)ifelse(is.infinite(q)|is.nan(q),0,q)),parent=den,count=num)
#		}
#		}))
#ts<-factor(tsort(x))
#d<-d[match(ts,factor(as.character(d[,1]))),]
# d<-data.frame(d)
#d[,2]<-as.numeric(as.character(d[,2]))
#return(d)
#})
#.export2LabKey<-function(obj,popnames,outpath){
#	stats<-getPopStats(obj)
#	rownames(stats)<-popnames;
#	popnames.t<-sapply(popnames,function(x)paste(strwrap(x,width=getOption("width")*0.6),collapse="\n"))
#	tsvout<-data.frame()
#	for(i in 2:length(popnames)){
#		plotfile<-paste(paste(file.path(outpath),basename(uuid()),sep="/"),"png",sep=".")
#		png(plotfile);
#		selectMethod("plot",signature=c("flowClustTree","character"))(obj,tsort(obj)[i]);
#		title(popnames.t[i]);
#		this.pop<-stats[i,1]
#		model<-getModel(obj,tsort(obj)[i])
#		xdim<-NA
#		ydim<-NA
#		if(!is.null(model)){
#			mu<-try(model@mu[which(attr(model@prior,"popnames")%in%this.pop),,drop=FALSE])
#			points(mu,pch="X");
#			#text(x=mu[1],y=mu[2],signif(stats[i,2],4))
#			xdim<-model@varNames[1]
#			ydim<-model@varNames[2]
#		}
#		dev.off();
#		tsvout<-rbind(tsvout,
#		rbind(c(SampleName=identifier(getData(obj,nodes(obj)[1])),PopName=popnames[i],Statistic="model.proportion",Value=stats[i,]$model.proportion.events,PlotFile=plotfile,xdim=xdim,ydim=ydim),
#		c(SampleName=identifier(getData(obj,nodes(obj)[1])),PopName=popnames[i],Statistic="flowCore.proportion",Value=stats[i,]$proportion.events,PlotFile=plotfile,xdim=xdim,ydim=ydim),
#		c(SampleName=identifier(getData(obj,nodes(obj)[1])),PopName=popnames[i],Statistic="parent.count",Value=stats[i,]$parent.events,PlotFile=plotfile,xdim=xdim,ydim=ydim),
#		c(SampleName=identifier(getData(obj,nodes(obj)[1])),PopName=popnames[i],Statistic="pop.count",Value=stats[i,]$count.events,PlotFile=plotfile,xdim=xdim,ydim=ydim)))
#	}
#	tsvout;
#}
#getPopNames<-function(x){
#	root<-nodes(x)[1];
#	d<-getPopStats(x)
#	sapply(as.character(d[,1]),function(x)strsplit(x,"\\.")[[1]][2])
#	newnms<-gsub("\\)","\\]",gsub("\\(","\\[",sapply(1:length(tsort(x)),function(i)paste(unlist(lapply(strsplit(RBGL::sp.between(x,root,tsort(x)[i])[[1]]$path_detail,"\\."),function(x)ifelse(length(x)==1,x[1],x[2]))),collapse="/"))))
#	newnms;
#}

#setMethod("plot",c("flowClustTree","GatingHierarchy"),function(x,y,node,level=0.9,pch=20,add=FALSE,...){
#	if(missing(node)){
#		stop("You need to provide an argument for node into the prior tree");
#	}
#	if(!is.character(node)){
#		stop("node must be of type character")
#	}
#		i<-match(node,bfs(x));
#		dims<-getDimensions(y,bfs(y@tree)[i])
#		d<-getData(x,node,parent=T)
#		dims<-dims[na.omit(match(parameters(d)@data$name,dims))]
#		if(!add){
#		plotGate(y,bfs(y@tree)[i],add=F)
#			selectMethod("plot",c("flowClustTree","character"))(x,bfs(x)[i],npoints=2000,add=T,level=level,pch=pch,subset=dims,main=bfs(x)[i],...)
#		}
#		if(class(getGate(y,bfs(y@tree)[i]))=="BooleanGate"){
#			plotGate(y,bfs(y@tree)[i],add=T,lwd=2,pch=20,...)
#		}else{
#			plotGate(y,bfs(y@tree)[i],add=T,lwd=2,...)
#		}
#})

# ===========================================================
# = Plot a gate or gates in nodes of a flowClust prior tree =
# = Note: you must call plot first!							=
# ===========================================================
#setMethod("plotGate",signature=c("flowClustPriorTree","character"),function(x,y,rev=F){
#	gates<-nodeData(x,y,"gate");
#	lapply(gates,function(x)lapply(x,function(x)polygon({nc<-ncol(x@boundaries);if(rev)x@boundaries[,nc:1] else x@boundaries},border="gray",lwd=1)))
#	
#})
# ============================
# = Plot flowClustTree graph =
# ============================
#setMethod("plot",c("flowClustTree","missing"),function(x,y,layoutType="dot",...){
#	renderGraph(Rgraphviz::layoutGraph(x,layoutType=layoutType,attrs=list(graph=list(rankdir="LR",page=c(8.5,11)),node=list(fixedsize=FALSE,fontsize=14,shape="rectangle"))))
#})
# =====================================================
# = S4 to fetch the data for the parent or given node =
# = in a flowClustTree object						  =
# =====================================================
#setMethod("getData",c("flowClustTree"),function(obj,y,parent=TRUE,...){
#	if(parent){
#		p<-.graph_parent(obj,y)
#		if(length(p)==0){
#			message("Node ", y," has no parent. Returning self.");
#			return(nodeData(obj,y,"population")[[1]])
#		}
#		return(nodeData(obj,p,"population")[[1]])
#	}else{
#		return(nodeData(obj,y,"population")[[1]])
#	}
#})
# =============================================
# = Fetch the model from a flowClustTree node =
# =============================================
#setGeneric("getModel",function(tree,name){
#	standardGeneric("getModel");
#})
# ==============================================
# = Fetcth the model from a flowClustTree node =
# ==============================================
#setMethod("getModel",signature("flowClustTree","character"),function(tree,name){
#	e<-graph::nodeData(tree,name,"model")[[1]]
#	if(is.environment(e)){
#		return(get("model",envir=e));
#	}else{
#		return(NULL)
#	}
#})

# =========================================================================================================
# = Map function for flowClust model with bayes prior. Assigns events to clusters based on quantile. =
# = Assigns outliers conditional on cluster membership, not globally. Important for correclty assigning   =
# = Rare cell populations																				  =
# =========================================================================================================
.fcbMap<-function(x,quantile){
	## TODO This needs to make use of the prior nu0.
	p<-length(x@varNames)
	qq<-qf(quantile,p,x@nu)
	qq<-(x@nu+p)/(x@nu+p*qq)
	(apply(x@z*apply(x@u,2,function(y)ifelse((y<qq),0,1)),1,function(x)ifelse(all(x==0),0,which.max(x))))
}

# =======================================================
# = fit a flowClust model to a node in a tree of priors.=
# = and a gating set of data							=
# = Input: GatingSet, prior tree, node index in			=
# = topological sort order.	The node index should be	=
# = the PARENT of the nodes whose gates you wish to fit =
# =======================================================
#flowClustBayesTree<-function(gatingset,priortree,nodeindex,...){
#	if(length(getChildren(gatingset[[1]],getNodes(gatingset[[1]],tsort=T)[nodeindex]))==0){
#		stop("node index ",nodeindex," is out of bounds. Node ",getNodes(gatingset[[1]],tsort=T)[nodeindex]," has no children");
#	}
#	if(class(gatingset)!="GatingSet"){
#		stop("gatingset must be of class GatingSet but got class ",class(gatingset));
#	}
#	if(class(priortree)!="flowClustPriorTree"){
#		stop("gatingset must be of class flowClustPriorTree but got class ",class(gatingset));
#	}
#	if(!is.numeric(nodeindex)){
#		stop("nodeindex should be numeric but got class ",class(nodeindex));
#	}
#	res<-sapply(1:length(gatingset),function(i)
#			flowClustBayes(x=getData(gatingset[[i]],nodeindex,tsort=TRUE),
#				prior=flowClust:::.combinePriorsForChildNodes(priortree,RBGL::tsort(priortree)[nodeindex]),...)
#			)
#	return(res);
#}

# =======================================================================================
# = Fetch the indices of cluser membership for a flowClustPriorTree, given a node name. =
# =======================================================================================
#setMethod("getIndices",signature=c("flowClustTree","character"),function(obj,y){
#	.getIndices(obj,y);
#})
 # =================================================================================================
 # = Code to get the full length class membership indices for flowClustTree used in F-measure code =
 # =================================================================================================
#.getIndices<-function(obj,y){
#	if(class(nodeData(obj,y,"gate")[[1]])=="BooleanGate"){
#		indices<-.getBooleanGateIndices(obj,y);
#		return(indices);
#	}else{
#		root<-nodes(obj)[1]
#		path<-RBGL::sp.between(obj,root,y)[[1]]$path_detail;
#		indices<-rep(TRUE,nrow(getData(obj,root)))
#		for(y in path[2:length(path)]){
#			if(class(nodeData(obj,y,"gate")[[1]])!="BooleanGate"){
#				m<-getModel(obj,y);
#				indices[which(indices)]<-indices[which(indices)]&(Map(m)==match(y,attr(m@prior,"popnames")))
#			}else{
#				indices<-indices&.getBooleanGateIndices(obj,y);
#			}
#		}
#	}
#	#indices[which(is.na(indices))]<-FALSE
#	return(indices);
#}
# ===========================================================
# = Code to get the indices for a boolean gate, full length =
# ===========================================================
#.getBooleanGateIndices<-function(obj,y){
#	#common ancestor
#	bg<-nodeData(obj,y,"gate")[[1]]
#	sps<-sapply(bg$ref,function(x){
#		RBGL::sp.between(obj,nodes(obj)[1],x)
#	})
#	#ca.ind<-min(unlist(lapply(sps,function(x)x$length+1)))
#	#ca<-sps[[1]]$path_detail[min(unlist(lapply(sps,function(x)x$length+1)))]
#	ca<-tsort(obj)[2]
#	ca.ind<-2;
#	amodel<-getModel(obj,ca)
#	indices<-list();
#		for(j in 1:length(sps)){
#			indices[[j]]<-rep(TRUE,dim(amodel@z)[1])
#			#indices[[j]]<-indices[[j]][which(indices[[j]]==T)]
#			for(k in (ca.ind):(sps[[j]]$length+1)){ 
#				this<-sps[[j]]$path_detail[k]
#				ti<-match(this,attr(getModel(obj,this)@prior,"popnames"));
#				indices[[j]][which(indices[[j]]==TRUE)]<-indices[[j]][which(indices[[j]]==TRUE)]&Map(getModel(obj,this))==ti
#			}	
#		}
#		query<-vector("list",length(indices));
#		query<-lapply(query,function(x)character())
#		for(j in 1:length(bg$ref)){
#			query[[j]]<-paste("indices[[",j,"]]",sep="")
#			query[[j]]<-paste(bg$v[j],query[[j]],sep="")
#		}
#		for(j in 1:length(bg$v2)){
#			query[[j+1]]<-(paste(unlist(query)[j:(j+1)],collapse=bg$v2[j]))
#		}
#		query<-query[[length(query)]];
#		#These are the indices of events included in the boolean gate relative to the parent of the component gates.
#		indices<-eval(parse(text=query))
#		indices[is.na(indices)]<-FALSE
#		return(indices);
#}
# ======================
# = Get the leaf nodes =
# ======================
#.getLeafNodes<-function(obj){
#	if(class(obj)=="GatingHierarchy")
#		obj<-obj@tree
#	names(which(unlist(lapply(adj(obj,tsort(obj)),function(x)length(x)==0))))
#}
#getCLROutput<-function(x){
#	inds<-sapply(.getLeafNodes(x),function(y)getIndices(x,y))
#	inds[is.na(inds)]<-FALSE
#	inds;
#}
# =======================================================
# = flowClust Model Fitting Wrapper for use with Priors =
# = Don't need to pass varNames. They will be taken     =
# = from the prior specification, same for K and 		=
# = usePrior.
# =======================================================
# flowClustBayes<-function(x, expName = "Flow Experiment", 
# 		varNames = if(class(prior)=="flowClustPrior") colnames(prior$Mu0) else NULL, K, 
# 	    B = 1000, tol = 1e-05, nu = 4, lambda = 1, nu.est = 0, trans = 1, 
# 	    min.count = 10, max.count = 10, min = NULL, max = NULL, level = 0.9, 
# 	    u.cutoff = NULL, z.cutoff = 0, randomStart = 10, B.init = B, 
# 	    tol.init = 0.01, seed = 1, criterion = "BIC", control = NULL, 
# 	    prior = NULL,usePrior="yes",w0=c(5,1,0.5,0.01),nu0=NULL,rare.thresh=c(0.02,0.15,0.35),pi0=NULL,...){
# 			K<-nrow(prior[[1]]$Mu0);
# 			if(class(prior)=="flowClustPriorList"){
# 				results<-vector("list",length(prior))
# 				#Loop through all the priors and fit the models.
# 				for(i in 1:length(prior)){
# 					varNames<-colnames(prior[[i]]$Mu0)
# 					#pi0 is prior # of observations for the data
# 					if(!is.null(pi0)){
# 						#if pi0 is null, then use the empirical bayes estimate of the prior cluster weights.
# 						#otherwise rescale the prior$w0, (the alpha hyperparameters)
# 						prior[[i]]$w0<-prior[[i]]$w0*pi0/sum(prior[[i]]$w0)
# 					}
# 					#nu0.hat and w0.hat are set to some defaults.
# 					nu0.hat<-prior[[i]]$nu0
# 					w0.hat<-rep(w0[4],K);
# 					
# 					#for each prior/model
# 					for(k in 1:nrow(prior[[i]]$Mu0)){
# 							#If the weights are less than the first threshold
# 							if(prior[[i]]$w0[k]/sum(prior[[i]]$w0)<rare.thresh[1]){
# 								#and if nu0 (the weighting factor for Lambda0) is provided
# 								if(!is.na(nu0[1])){
# 									#set the weighting factor for lambda0 (nug0.star in the paper) for the k'th prior to the given value.
# 						 			nu0.hat[k]<-nu0[1]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)
# 								}
# 								#check the second interval
# 							}else if(prior[[i]]$w0[k]/sum(prior[[i]]$w0)>rare.thresh[1]&prior[[i]]$w0[k]/sum(prior[[i]]$w0)<rare.thresh[2]){
# 								if(!is.na(nu0[2])){
# 									#set the lambda0 weight
# 									nu0.hat[k]<-nu0[2]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)
# 								}
# 							}
# 							#Finally, reweight Lambda0	
# 								prior[[i]]$Lambda0[k,,]<-((nu0.hat[k]-ncol(prior[[i]]$Mu0)-1)*prior[[i]]$Lambda0[k,,])/(prior[[i]]$nu0[k]-ncol(prior[[i]]$Mu0)-1)
# 								#Don't forget to update nu0.
# 								prior[[i]]$nu0[k]<-nu0.hat[k]
# 							
# 							#check the intervals for the omega0 weight.	
# 							if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[1]){
# 								w0.hat[k]<-w0[1]
# 							}
# 							if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[2]&(prior[[i]]$w0[k]/sum(prior[[i]]$w0))>=rare.thresh[1]){
# 								w0.hat[k]<-w0[2]
# 							}
# 							if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))>=rare.thresh[2]&(prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[3]){
# 								w0.hat[k]<-w0[3]
# 							}
# 							#reweight omega0
# 								d<-ncol(prior[[i]]$Omega[k,,]);
# 								ne<-nrow(x); #total events
# 								Nk<-ne*prior[[i]]$w0[k]/sum(prior[[i]]$w0); #Ng based on prior cluster proportion
# 								w0.hat[k]<-w0.hat[k]*Nk*(det(prior[[i]]$Omega0[k,,])/det(prior[[i]]$Lambda0[k,,]/(prior[[i]]$nu0[k]-d-1)))^(1/d);
# 								prior[[i]]$Omega0[k,,]<-prior[[i]]$Omega0[k,,]*(1/(w0.hat[k]))
# 					}
# 					results[[i]]<-flowClust(x=x,expName=expName,varNames=varNames,K=nrow(prior[[i]]$Mu0),
# 					B=B,tol=tol,nu=nu,lambda=lambda,nu.est=nu.est,trans=trans,
# 					min.count=min.count,max.count=max.count,min=min,max=max,level=level,
# 					u.cutoff=u.cutoff,z.cutoff=z.cutoff,randomStart=randomStart,B.init=B.init,
# 					tol.init=tol.init,seed=seed,criterion=criterion,control=control,prior=prior[[i]],
# 					usePrior=usePrior,...)
# 				}
# 				return(results)
# 			}else{
# 				i<-1;
# 				if(!is.null(pi0)){
# 					prior[[i]]$w0<-prior[[i]]$w0*pi0/sum(prior[[i]]$w0)
# 				}
# 				nu0.hat<-prior[[i]]$nu0
# 				w0.hat<-rep(w0[4],K);
# 				
# 					for(k in 1:nrow(prior[[i]]$Mu0)){
# 								if(prior[[i]]$w0[k]/sum(prior[[i]]$w0)<rare.thresh[1]){
# 										nu0.hat[k]<-nu0[1]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)
# 								}else if(prior[[i]]$w0[k]/sum(prior[[i]]$w0)>rare.thresh[1]&prior[[i]]$w0[k]/sum(prior[[i]]$w0)<rare.thresh[2]){
# 										nu0.hat[k]<-nu0[2]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)
# 								}
# 								prior[[i]]$Lambda0[k,,]<-((nu0.hat[k]-ncol(prior[[i]]$Mu0)-1)*prior[[i]]$Lambda0[k,,])/(prior[[i]]$nu0[k]-ncol(prior[[i]]$Mu0)-1)
# 									prior[[i]]$nu0[k]<-nu0.hat[k]
# 								if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[1]){
# 									w0.hat[k]<-w0[1]
# 								}
# 								if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[2]&(prior[[i]]$w0[k]/sum(prior[[i]]$w0))>=rare.thresh[1]){
# 									w0.hat[k]<-w0[2]
# 								}
# 								if((prior[[i]]$w0[k]/sum(prior[[i]]$w0))>=rare.thresh[2]&(prior[[i]]$w0[k]/sum(prior[[i]]$w0))<rare.thresh[3]){
# 									w0.hat[k]<-w0[3]
# 								}
# 									d<-ncol(prior[[i]]$Omega[k,,]);
# 									ne<-nrow(x);
# 									Nk<-ne*prior[[i]]$w0[k]/sum(prior[[i]]$w0);
# 									w0.hat[k]<-w0.hat[k]*Nk*(det(prior[[i]]$Omega0[k,,])/det(prior[[i]]$Lambda0[k,,]/(prior[[i]]$nu0[k]-d-1)))^(1/d);
# 									prior[[i]]$Omega0[k,,]<-prior[[i]]$Omega0[k,,]*(1/(w0.hat[k]))
# 						}
# 				#Just return the model fitted using the prior
# 				
# 				return(flowClust(x=x,expName=expName,varNames=varNames,K=nrow(prior$Mu0),
# 				B=B,tol=tol,nu=nu,lambda=lambda,nu.est=nu.est,trans=trans,
# 				min.count=min.count,max.count=max.count,min=min,max=max,level=level,
# 				u.cutoff=u.cutoff,z.cutoff=z.cutoff,randomStart=randomStart,B.init=B.init,
# 				tol.init=tol.init,seed=seed,criterion=criterion,control=control,prior=prior,
# 				usePrior=usePrior,...));
# 			}
# 				
# }


##Rewrite flowClustBayes to use a variable length interval list
# =======================================================
# = flowClust Model Fitting Wrapper for use with Priors =
# = Don't need to pass varNames. They will be taken     =
# = from the prior specification, same for K and 		=
# = usePrior.
# =======================================================
# TODO Try setting  a default minimum number of observations for Lambda0, Omega0
# TODO Try setting the equivalent observations for each cluster to the same value.
#flowClustBayes<-function(x, expName = "Flow Experiment", 
#		varNames = if(class(prior)=="flowClustPrior") colnames(prior$Mu0) else NULL, K, 
#	    B = 1000, tol = 1e-05, nu = 4, lambda = 1, nu.est = 0, trans = 1, 
#	    min.count = 10, max.count = 10, min = NULL, max = NULL, level = 0.9, 
#	    u.cutoff = NULL, z.cutoff = 0, randomStart = 10, B.init = B, 
#	    tol.init = 0.01, seed = 1, criterion = "BIC", control = NULL, 
#	    prior = NULL,usePrior="yes",w0=c(5,1,0.5,0.01),nu0=NULL,rare.thresh=c(0.02,0.15,0.35),pi0=NULL,...){
#			if(class(prior)=="flowClustPriorList"){
#				results<-vector("list",length(prior))
#				#Loop through all the priors and fit the models.
#				for(i in 1:length(prior)){
#					K<-nrow(prior[[i]]$Mu0);
#					message("Fitting ",paste(attr(prior[[i]],"popnames"),collapse=" "))
#					varNames<-colnames(prior[[i]]$Mu0)
#					#pi0 is prior # of observations for the data
#					if(!is.null(pi0)){
#						#if pi0 is null, then use the empirical bayes estimate of the prior cluster weights.
#						#otherwise rescale the prior$w0, (the alpha hyperparameters)
#						prior[[i]]$w0<-prior[[i]]$w0*pi0/sum(prior[[i]]$w0)
#					}
#					#if(attr(prior[[i]],"popnames")[1]=="45.Activated CD4 T cells [CD25]"){
#					#	browser();
#					#}
#					#nu0.hat and w0.hat are set to some defaults.
#					nu0.hat<-prior[[i]]$nu0
#					w0.hat<-rep(0.01,K); #primarily data-driven
#					#for each component
#					for(k in 1:nrow(prior[[i]]$Mu0)){
#						
#							#Scan the intervals
#							P<-prior[[i]]$w0[k]/sum(prior[[i]]$w0)
#							wh<-rev(c(ifelse(any(P<=rare.thresh),min(which(P<=rare.thresh)),length(rare.thresh)+1),ifelse(any(P>rare.thresh),max(which(P>rare.thresh)),0)))[2]
#							#nu0 is a scaling factor for LambdaG0
#							if(!is.null(nu0)){
#							if(wh<length(nu0)){
#								nu0.hat[k]<-nu0[wh]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)								
#							}
#						}
#							#Ensure that the weight doesn't make Lambda0 negative.
#							nu0.hat[k]<-max(prior[[i]]$nu0[k],nu0.hat[k])
#							#Reweight Lambda0
#							prior[[i]]$Lambda0[k,,]<-((nu0.hat[k]-ncol(prior[[i]]$Mu0)-1)*prior[[i]]$Lambda0[k,,])/(prior[[i]]$nu0[k]-ncol(prior[[i]]$Mu0)-1)
#							#Don't forget to update nu0.
#							prior[[i]]$nu0[k]<-nu0.hat[k]
#							if(!is.null(w0)){
#							if(wh<=length(w0)){
#								w0.hat[k]<-w0[wh]
#							}}
#							#Reweight Omega0
#							if(class(prior[[i]]$Omega[k,,])!="matrix"){
#								d<-length(prior[[i]]$Omega[k,,])
#							}else{
#								d<-ncol(prior[[i]]$Omega[k,,]);
#							}
#							ne<-nrow(x); #total events
#							Nk<-ne*prior[[i]]$w0[k]/sum(prior[[i]]$w0); #Ng based on prior cluster proportion
#							w0.hat[k]<-w0.hat[k]*Nk*(.det(prior[[i]]$Omega0[k,,])/.det(prior[[i]]$Lambda0[k,,]/(prior[[i]]$nu0[k]-d-1)))^(1/d);
#							prior[[i]]$Omega0[k,,]<-prior[[i]]$Omega0[k,,]*(1/(w0.hat[k]))
#							
#					}
#					#if(attr(prior[[i]],"popnames")[1]=="45.Activated CD4 T cells [CD25]"){
#					#	debug(.flowClustK)
#					#}
#					results[[i]]<-flowClust(x=x,expName=expName,varNames=varNames,K=nrow(prior[[i]]$Mu0),
#					B=B,tol=tol,nu=nu,lambda=lambda,nu.est=nu.est,trans=trans,
#					min.count=min.count,max.count=max.count,min=min,max=max,level=level,
#					u.cutoff=u.cutoff,z.cutoff=z.cutoff,randomStart=randomStart,B.init=B.init,
#					tol.init=tol.init,seed=seed,criterion=criterion,control=control,prior=prior[[i]],
#					usePrior=usePrior,...)
#				}
#				return(results)
#			}else{
#				i<-1;
#				if(!is.null(pi0)){
#					prior[[i]]$w0<-prior[[i]]$w0*pi0/sum(prior[[i]]$w0)
#				}
#				message("Fitting ",paste(attr(prior,"popnames"),collapse=" "))
#				#set some defaults
#				nu0.hat<-prior[[i]]$nu0
#				w0.hat<-rep(0.01,K); #data driven
#				
#					for(k in 1:nrow(prior[[i]]$Mu0)){
#						#Scan the intervals
#						P<-prior[[i]]$w0[k]/sum(prior[[i]]$w0)
#						wh<-rev(c(ifelse(any(P<=rare.thresh),min(which(P<=rare.thresh)),length(rare.thresh)+1),ifelse(any(P>rare.thresh),max(which(P>rare.thresh)),0)))[2]
#						if(!is.null(nu0)){
#						if(wh<length(nu0)){
#							nu0.hat[k]<-nu0[wh]*prior[[i]]$w0[k]/sum(prior[[i]]$w0)*nrow(x)	#nu0 is a scaling factor for the number of observations.							
#						}
#					}
#						#Ensure that the weight doesn't make Lambda0 negative.
#						nu0.hat[k]<-max(prior[[i]]$nu0[k],nu0.hat[k])
#						#Reweight Lambda0
#						prior[[i]]$Lambda0[k,,]<-((nu0.hat[k]-ncol(prior[[i]]$Mu0)-1)*prior[[i]]$Lambda0[k,,])/(prior[[i]]$nu0[k]-ncol(prior[[i]]$Mu0)-1)
#						#Don't forget to update nu0.
#						prior[[i]]$nu0[k]<-nu0.hat[k]
#						if(!is.null(w0)){
#						if(wh<=length(w0)){
#							w0.hat[k]<-w0[wh]
#						}
#						}
#						#Reweight Omega0
#						if(class(prior[[i]]$Omega[k,,])!="matrix"){
#							d<-length(prior[[i]]$Omega[k,,])
#						}else{
#							d<-ncol(prior[[i]]$Omega[k,,]);
#						}
#						ne<-nrow(x); #total events
#						Nk<-ne*prior[[i]]$w0[k]/sum(prior[[i]]$w0); #Ng based on prior cluster proportion
#						w0.hat[k]<-w0.hat[k]*Nk*(.det(prior[[i]]$Omega0[k,,])/.det(prior[[i]]$Lambda0[k,,]/(prior[[i]]$nu0[k]-d-1)))^(1/d);
#						prior[[i]]$Omega0[k,,]<-prior[[i]]$Omega0[k,,]*(1/(w0.hat[k]))
#						
#						}
#				#Just return the model fitted using the prior
#				#if(attr(prior,"popnames")[1]=="45.Activated CD4 T cells [CD25]"){
#				#	debug(.flowClustK)
#				#}
#				return(flowClust(x=x,expName=expName,varNames=varNames,K=nrow(prior$Mu0),
#				B=B,tol=tol,nu=nu,lambda=lambda,nu.est=nu.est,trans=trans,
#				min.count=min.count,max.count=max.count,min=min,max=max,level=level,
#				u.cutoff=u.cutoff,z.cutoff=z.cutoff,randomStart=randomStart,B.init=B.init,
#				tol.init=tol.init,seed=seed,criterion=criterion,control=control,prior=prior,
#				usePrior=usePrior,...));
#			}
#				
#}
##Need a wrapper for "det" to handle 1D numeric 1x1 matrix and return the correct result
.det<-function(x){
	if(class(x)=="numeric"){
		x<-as.matrix(x);
	}
	as.matrix(det(x));
}
# ===================================================================================================
# = Compare population statistics between a flowClustTree fit to a GatingSet and the GatingSet      =
# ===================================================================================================
#comparePopStats<-function(x,y){
#	if(!(all(unlist(lapply(x,function(x)inherits(x,"flowClustTree"))))&inherits(y,"GatingSet")&length(x)==length(y))){
#		stop("comparePopStats: invalid arguments. provide a list of flowClustTree and a GatingSet, with samples in the same order");
#	}
#	a<-do.call(cbind,lapply(x,function(x)getPopStats(x)[,2]))
#	b<-do.call(cbind,lapply(y,function(x)getPopStats(x)[,2]))
#	n<-getPopStats(x[[1]])[,1]
#	d<-data.frame(cv.x=apply(a,1,function(x)sd(x)/mean(x)),cv.y=apply(b,1,function(x)sd(x)/mean(x)),diff=rowMeans(a)-rowMeans(b))
#	rownames(d)<-n;
#	d
#}

