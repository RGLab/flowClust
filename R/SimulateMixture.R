SimulateMixture<-function(N,nu=100,mu,sigma,w)
{
	count<-0
	# Number of clusters
	K<-nrow(mu)
	# dimension
	p<-ncol(mu)
	y<-matrix(0,N,p)

	for(i in 1:N)
	{
		u<-runif(1)
		SumW<-w[1]
		count<-1
		while((count<K) & (u>SumW))
		{
			count<-count+1
			SumW<-SumW+w[count]			
		}
		y[i,]<-rmt(n = 1, mean=mu[count,], sigma[count,,], df=nu)
	}
	y
}
