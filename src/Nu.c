#include "flowClust.h"

double NuGradient(double x, void *params)
{
    struct Nu_params *p = (struct Nu_params *) params;     
    gsl_vector *SumZ = p->SumZ;
    gsl_vector *SumZU = p->SumZU;    
    gsl_vector *SumZlogU = p->SumZlogU;
    int K1 = p->K;

    int k=0, K=SumZ->size;
    int K0=0;
    double Tmp=0;
    double logLike=0;   // derivative of loglikelihood wrt nu
    // Rprintf("K=%d  nu=%lf\n",K1,x);

    Tmp=gsl_sf_log(x/2)-gsl_sf_psi(x/2)+1;
    if(K1>-1) K0=K1;
    for(k=K0;k<K;k++)
    {
        logLike+=gsl_vector_get(SumZlogU,k)-gsl_vector_get(SumZU,k);
        logLike+=gsl_vector_get(SumZ,k) * Tmp;
        if(K1>-1) break;
    }
    
    // Rprintf("logLike = %lf\n",logLike);
    return(logLike);
}


double NuLikelihood(double x, void *params)
{
    struct Nu_ECME_params *p = (struct Nu_ECME_params *) params;     
    gsl_matrix *YTrans = p->YTrans;
    gsl_matrix *Mu = p->Mu;
    gsl_matrix *Precision = p->Precision;
    gsl_matrix *logY = p->logY;
    gsl_vector *W = p->W;
    double *lambda = p->lambda;
    int *transform = p->transform;

    int i=0, ly=YTrans->size1;    
    int j=0, py=Mu->size2;
    int k=0, K=Mu->size1;
	gsl_vector_view rowYTrans, subRowYTrans, rowMu, rowPrecision;
	gsl_matrix_view matrixPrecision;
    double logJacobian=0;
    double like=0;    // Sum of (W * density function) wrt one observation
    double logLike=0;   // actual-data loglikelihood
    // Rprintf("nu=%lf\n",x);

	for(i=0;i<ly;i++)
	{        
        if(*transform==1 || (*transform==0 && *lambda!=1)) 
        {
    		for(j=0;j<py;j++)
		    	logJacobian+=(*lambda-1)*gsl_matrix_get(logY,i,j);
		}
		rowYTrans=gsl_matrix_row(YTrans,i);
        like=0;
		for(k=0;k<K;k++)
		{
			rowMu=gsl_matrix_row(Mu,k);
			rowPrecision=gsl_matrix_row(Precision,k);
			matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,py,py);
			/* Here I assume the cholesky decomposition has been done */
			/* Compute W * density function */
			if(*transform<=1)
    			like+=gsl_vector_get(W,k) * gsl_ran_mvnt_pdf(&rowYTrans.vector,&rowMu.vector,&matrixPrecision.matrix,x,1,0);
    	    else
    	    {
                logJacobian=0;
                for(j=0;j<py;j++)
                    logJacobian+=gsl_matrix_get(logY,i,j);
                logJacobian*=lambda[k]-1;
    	    
                subRowYTrans=gsl_vector_subvector(&rowYTrans.vector,k*py,py);
                like+=gsl_vector_get(W,k) * gsl_ran_mvnt_pdf(&subRowYTrans.vector,&rowMu.vector,&matrixPrecision.matrix,x,1,0) * gsl_sf_exp(logJacobian);            
    	    }
        }
		logLike+=log(like);
    }
	if(*transform==1 || (*transform==0 && *lambda!=1)) 
		logLike+=logJacobian;
    
    // Rprintf("logLike = %lf\n",logLike);
    return(-logLike);
}
