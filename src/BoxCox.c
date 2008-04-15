#include "flowClust.h"

double BoxCoxGradient(double x, void *params)
{
    struct BoxCox_params *p = (struct BoxCox_params *) params;     
    gsl_matrix *Mu  = p->Mu;
    gsl_matrix *Precision  = p->Precision;
    gsl_vector *W  = p->W;
    gsl_matrix *Weight = p->Weight;
    gsl_matrix *Z = p->Z;
    gsl_matrix *Y = p->Y;
    gsl_matrix *WeightedY = p->WeightedY;    
    gsl_matrix *DiagOne = p->DiagOne;
    gsl_vector *SumWZ = p->SumWZ;                
    gsl_vector *SumZ = p->SumZ;    

    int k=0,K=Mu->size1;
    int i=0,ly=Y->size1;    
    int j=0,py=Mu->size2;
    /*    printf("lambda=%lf\n",x);*/
    
    gsl_vector_view Ylong, Y1long, row, col;
    gsl_matrix_view YY, YY1, Prow;
    double logLike=0,Yoriginal=0,logJacobian=0,Tmp=0;
    gsl_vector *SMu = gsl_vector_alloc(py), *YTmp = gsl_vector_alloc(ly);
    gsl_matrix *Yo = gsl_matrix_alloc(ly,py);
    gsl_matrix *WeightedY1 = gsl_matrix_alloc(ly,K*py);

    /* Change all the data to include the transformation */
    for(i=0;i<ly;i++)
    {
        for(j=0;j<py;j++)
        {
            /* Retransform Y */
            Yoriginal=gsl_matrix_get(Y,i,j);
            gsl_matrix_set(Yo,i,j,(sgn(Yoriginal)*pow(fabs(Yoriginal),x)-1.0)/x);
            for(k=0;k<K;k++)
            {
                gsl_matrix_set(WeightedY,i,k*py+j,gsl_matrix_get(Weight,i,k)* gsl_matrix_get(Yo,i,j));
                gsl_matrix_set(WeightedY1,i,k*py+j,gsl_matrix_get(Weight,i,k)* ( sgn(Yoriginal)*pow(fabs(Yoriginal),x) * (log(fabs(Yoriginal))*x-1) +1) /(x*x));
            }
            logJacobian+=log(fabs(Yoriginal));
        }
    }
    
    /* Compute the estimate of Mu */  
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Z, Yo, 0, Mu);

    for(k=0;k<K;k++)
    { 
        row=gsl_matrix_row(Precision,k);
        Prow=gsl_matrix_view_vector(&row.vector,py, py);            
        /* Update Mu*/
        row=gsl_matrix_row(Mu,k);
        gsl_blas_dscal(1./gsl_vector_get(SumWZ,k),&row.vector);
		/* Update Precision (cluster specific) */
        /* Compute the estimate of Sigma inverse */                         
        YY=gsl_matrix_submatrix(WeightedY, 0, k*py, ly, py);  
        up_date_precision(&YY.matrix, &row.vector, &Prow.matrix, /*&col.vector,*/ gsl_vector_get(SumZ,k), gsl_vector_get(SumWZ,k), DiagOne);
        /* Update Mixing Proportions */
        gsl_vector_set(W,k,gsl_vector_get(SumZ,k)/ly);

        YY1=gsl_matrix_submatrix(WeightedY1, 0, k*py, ly, py);               
        gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &Prow.matrix, &YY1.matrix);        
        gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &Prow.matrix, &YY.matrix);                
		for(i=0;i<ly;i++)
		{
			Ylong=gsl_vector_view_array(WeightedY->data+i*py*K+k*py,py);
			Y1long=gsl_vector_view_array(WeightedY1->data+i*py*K+k*py,py);
			gsl_blas_ddot(&Ylong.vector, &Y1long.vector, &Tmp);
			logLike-=Tmp;
		}
	
		gsl_vector_memcpy(SMu, &row.vector);
        gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &Prow.matrix, SMu);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &YY1.matrix, SMu, 0.0, YTmp);
        col=gsl_matrix_column(Weight,k);            
        gsl_blas_ddot(YTmp, &col.vector, &Tmp);
		logLike+=Tmp;               
    }
    
    logLike+=logJacobian;     
	gsl_vector_free(SMu);
    gsl_vector_free(YTmp);
    gsl_matrix_free(Yo);
    gsl_matrix_free(WeightedY1);
    
    /* If the LogLike is not finite, return the largest number possible */
	/*    if((gsl_isinf(logLike)!=0) || (gsl_isnan(logLike)==1))
    {
        logLike=-FLT_MAX;
    }
	*/    
	/*	printf("logLike=%lf\n",logLike);	*/    
	return(logLike);
}


int sgn(double x)
{
	if (x >= 0) return 1; else return -1;
}
