#include "flowClust.h"
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CMethodDef CEntries[] = {
  {"flowClust", (DL_FUNC)&flowClust, 20},
  {NULL, NULL, 0}
};

// intended to comment out R_init_flowClust. We do not encourage users to call the C routines in R.  We have already defined the corresponding R functions for users to use.
/* void R_init_flowClust(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    // R_useDynamicSymbols(dll, FALSE);
} */

void flowClust(double *y, int *ly, int *py, double *mu, double *precision, double *w, double *z, double *u, double *lambda, int *label, double *uncertainty, double *u_cutoff, double *z_cutoff, int *flagOutliers, int *K, double *nu, int *B, double *tol, int *transform, double *logLike)
{
	
	int i=0,j=0,k=0;
//	const gsl_rng_type * T;
//	gsl_rng * r;
	double logLikeOld=0,Diff=100.0;
	gsl_vector *SumZ=gsl_vector_calloc(*K), *SumWZ=gsl_vector_calloc(*K);
	gsl_matrix *DiagOne=gsl_matrix_calloc(*py,*py);
	gsl_matrix *YTrans=gsl_matrix_alloc(*ly,*py), *WeightedY=gsl_matrix_calloc(*ly,*K**py);
	gsl_matrix_view Y,Mu, Z, Precision, Prow, YY, Weight;
	gsl_vector_view W,row,row1,row3,row4,col,subRow;
	int iter=0;    
	const gsl_root_fsolver_type *S; 
	gsl_root_fsolver *s; 
	double x_lo= .1, x_hi=1, lambdaHat=0;
	double xLow=.1, xUp=1;
	int status=0, iterSolve=0, iterSolveMax=50;
	struct BoxCox_params params;
	gsl_function F;
	

	/** Turn off the error handler **/
	gsl_set_error_handler_off();

	/* This is actually not need but I leave it for now if we want to look at stochastic EM in the future
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(T); */

	S=gsl_root_fsolver_brent;
	s=gsl_root_fsolver_alloc(S);
	F.function=&BoxCoxGradient;
	F.params=&params;

	/* Create matrix and vector views*/  
	Y=gsl_matrix_view_array(y,*ly,*py);
	/* Add a constant to all values to avoid problems with zeros */	
  gsl_matrix_add_constant(&Y.matrix, 0.000001111);
    /* Set xLow to -1 if all observations are positive */
//    if(gsl_matrix_min(&Y.matrix)>0) xLow=-1;
	Mu=gsl_matrix_view_array(mu,*K,*py);
	Z=gsl_matrix_view_array(z,*ly,*K);
	Weight=gsl_matrix_view_array(u,*ly,*K);
	W=gsl_vector_view_array(w,*K);
	Precision=gsl_matrix_view_array(precision,*K,*py**py);

	/* Initialization */
	gsl_matrix_set_identity(DiagOne);	

	for(i=0;i<*ly;i++)
	{
		for(j=0;j<*py;j++)
		{
			gsl_matrix_set(YTrans,i,j,(sgn(gsl_matrix_get(&Y.matrix,i,j))*pow(fabs(gsl_matrix_get(&Y.matrix,i,j)),*lambda)-1)/(*lambda));
		}
		row=gsl_matrix_row(WeightedY,i);
		row1=gsl_matrix_row(YTrans,i);        
		if(label[i]>0)
		{
			gsl_matrix_set(&Z.matrix,i,label[i]-1,1.0);
			/* I initialize the weights only for observations that have Z_i==1 */
			gsl_matrix_set(&Weight.matrix,i,label[i]-1,1.0);
			gsl_vector_set(SumZ,label[i]-1,gsl_vector_get(SumZ,label[i]-1)+1.0);
			gsl_vector_set(SumWZ,label[i]-1,gsl_vector_get(SumWZ,label[i]-1)+1.0);
		}
		for(k=0;k<*K;k++)
		{
			/* Only extract the weighted Y's for the kth cluster and reweight with new weights */
			subRow=gsl_vector_subvector(&row.vector, k**py, *py);
			gsl_vector_memcpy(&subRow.vector,&row1.vector);
			gsl_blas_dscal(gsl_matrix_get(&Weight.matrix,i,k), &subRow.vector);
		}
	}

	/* Initialization */
	/* Mu using BLAS */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans, 0, &Mu.matrix);
	for(k=0;k<*K;k++)
	{
		/* This could be used to do a cluster specific update of Mu*/
		/*  up_date_mu(&Y.matrix, &row.vector, &col.vector, gsl_vector_get(SumWZ,k));  */
		row=gsl_matrix_row(&Precision.matrix,k);
		Prow=gsl_matrix_view_vector(&row.vector,*py, *py);            
		row=gsl_matrix_row(&Mu.matrix,k);
		gsl_blas_dscal(1./gsl_vector_get(SumWZ,k),&row.vector);
		/* Update Precision (cluster specific) */
		YY=gsl_matrix_submatrix(WeightedY, 0, k**py, *ly, *py);      
/**/		col=gsl_matrix_column(&Weight.matrix,k);            
		up_date_precision(&YY.matrix, &row.vector, &Prow.matrix, /* &col.vector, */ gsl_vector_get(SumZ,k), gsl_vector_get(SumWZ,k), DiagOne);            
        /* Update Mixing Proportions */
		gsl_vector_set(&W.vector,k,gsl_vector_get(SumZ,k)/(*ly));
	}
	/* Initialize the log likelihood */
	*logLike=-FLT_MAX;

	/* EM algorithm */
	while((Diff>*tol) & (iter<*B))
	{		
		logLikeOld=*logLike;
		/* E step*/
		up_date_z_weight(&Y.matrix, YTrans, /*WeightedY,*/ &Mu.matrix, &Precision.matrix, &W.vector, &Z.matrix, &Weight.matrix, SumZ, SumWZ, *nu, *lambda, logLike, 0, *transform);      
		/* M step*/        
		/** Solve for the Box-Cox parameter **/
		params.Mu=&Mu.matrix;
		params.Precision=&Precision.matrix;
		params.W=&W.vector;        
		params.Weight=&Weight.matrix;
		params.Z=&Z.matrix;
		params.Y=&Y.matrix;
		params.WeightedY=WeightedY;
		params.SumWZ=SumWZ;
		params.SumZ=SumZ;
		params.DiagOne=DiagOne;        
		params.lambda=*lambda;

		F.params=&params;

		/* Only estimate lambda in the first 100 iterations */
		if(iter>=100)
			*transform=0;			
			
		if(*transform==1)
		{
			status=gsl_root_fsolver_set(s, &F, xLow, xUp);
/*			if(status!=0)
			{
				warning("Brent's algorithm returned an error, error code %d\n", status);				
			}
*/                
		}
		if((*transform==1) && (status==0))
		{
			iterSolve=0;            
			do
			{
				iterSolve++;
				status=gsl_root_fsolver_iterate(s);
				lambdaHat=gsl_root_fsolver_root(s);
				x_lo=gsl_root_fsolver_x_lower(s);
				x_hi=gsl_root_fsolver_x_upper(s);
				status=gsl_root_test_interval(x_lo, x_hi,0.001, 0);
			}
			while((status == GSL_CONTINUE) && (iterSolve < iterSolveMax));
/*			if((status!=0) && (iterSolve<iterSolveMax))
			{
				warning("Warning in the inner loop, Brent's algorithm returned an error, error code %d\n", status);
			}
*/			/* Set the new value of the transformation */
			*lambda=lambdaHat;        

		}
		else if((*transform==0) || ((*transform==1) && (status!=0)))
		{
			/* Need to rescale with the proper weight */
			for(i=0;i<*ly;i++)
			{
				row1=gsl_matrix_row(YTrans,i);
				row3=gsl_matrix_row(WeightedY,i);                    
				for(k=0;k<*K;k++)
				{                        
					row4=gsl_vector_subvector(&row3.vector, k**py, *py);
					gsl_vector_memcpy(&row4.vector,&row1.vector);                        
					gsl_blas_dscal(gsl_matrix_get(&Weight.matrix,i,k), &row4.vector);                        
				}
			}

			/* Mu using BLAS */
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans, 0, &Mu.matrix);        
			for(k=0;k<*K;k++)
			{
				row=gsl_matrix_row(&Precision.matrix,k);
				Prow=gsl_matrix_view_vector(&row.vector,*py,*py);            
				row=gsl_matrix_row(&Mu.matrix,k);
/**/				col=gsl_matrix_column(&Weight.matrix,k);            
				gsl_blas_dscal(1./gsl_vector_get(SumWZ,k),&row.vector);    
				
		        /* Update Precision (cluster specific) */
				YY=gsl_matrix_submatrix(WeightedY, 0, k**py, *ly, *py);
				up_date_precision(&YY.matrix, &row.vector, &Prow.matrix, /* &col.vector, */ gsl_vector_get(SumZ,k), gsl_vector_get(SumWZ,k), DiagOne);                        
                /* Update Mixing Proportions */
				gsl_vector_set(&W.vector,k, gsl_vector_get(SumZ,k)/(*ly));            
			}

		}
		Diff=fabs(*logLike-logLikeOld)/fabs(logLikeOld);
		iter++; 
	}
	printf("The EM required %d iterations\n",iter);
	printf("The tolerance is %g\n",Diff);

	/* One more iteration to compute the final weight and u's */
	up_date_z_weight(&Y.matrix, YTrans, /*WeightedY,*/ &Mu.matrix, &Precision.matrix, &W.vector, &Z.matrix, &Weight.matrix, SumZ, SumWZ, *nu, *lambda, logLike, 1, *transform);      

    /* Output Precision to be the covariance matrix */
	for(k=0;k<*K;k++)
	{
	    gsl_matrix_set_identity(DiagOne);
		row=gsl_matrix_row(&Precision.matrix,k);
		Prow=gsl_matrix_view_vector(&row.vector,*py, *py);            
        gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &Prow.matrix, DiagOne);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, &Prow.matrix);
    }

	/** Compute the labels and uncertainty **/
	for(i=0;i<*ly;i++)
	{
		row=gsl_matrix_row(&Z.matrix,i);
		uncertainty[i]=1-gsl_vector_max(&row.vector);
		label[i]=gsl_vector_max_index(&row.vector)+1;
		if(gsl_matrix_get(&Weight.matrix,i,label[i]-1)<*u_cutoff || gsl_matrix_get(&Z.matrix,i,label[i]-1)<*z_cutoff)
			flagOutliers[i]=1;
        /*	else
			flagOutliers[i]=0; */
	}
	
	
	gsl_vector_free(SumZ);
	gsl_vector_free(SumWZ);
	gsl_matrix_free(DiagOne);
	gsl_matrix_free(WeightedY);
	gsl_matrix_free(YTrans);
	gsl_root_fsolver_free(s);
}

// void up_date_mu(gsl_matrix *Y, gsl_vector *Mu, gsl_vector *Weight, double SumWZ)
// {
	/* Compute W^TY*/
	/* Weight contains the Z_{ik}W_{ik} */
	/* gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Weight, Y, 0, Mu); */
//	gsl_blas_dgemv(CblasTrans, 1.0, Y, Weight, 0, Mu);

	/* Divide by \sum_iW_ikZ_ik */
//	gsl_blas_dscal(1./SumWZ, Mu);
// }

void up_date_precision(gsl_matrix *Y, gsl_vector *Mu, gsl_matrix *Precision, /* gsl_vector *Weight, */ double SumZ, double SumWZ, gsl_matrix *DiagOne)
{
	int status=0;
	gsl_blas_dsyrk(CblasLower, CblasTrans, 1./SumZ, Y, 0.0, Precision);
	gsl_blas_dsyr(CblasLower, -SumWZ/SumZ, Mu, Precision);    
	/* Compute the precision matrix and its cholesky decomposition */
	/* Set DiaOne to the identity*/
	gsl_matrix_set_identity(DiagOne);
	/* Compute the inverse (precision) and its cholesky decomposition once for all */	
	/* Compute the cholesky decomposition of the covariance matrix */
	status=gsl_linalg_cholesky_decomp(Precision);
	if(status!=0)
	{
		error("\n The covariance matrix is near singular! \n Try running the program with a different initial configuration or less clusters \n");		
	}
	/* Compute L'^{-1} */
	gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, Precision, DiagOne);
	/* Compute Sigma=L'^{-1}L^{-1} */
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, Precision);
	/* Compute the cholesky decomposition of the precision matrix */
	status=gsl_linalg_cholesky_decomp(Precision);
	if(status!=0)
	{
		error("\n The covariance matrix is near singular! \n Try running the program with a different initial configuration or less clusters \n");		
	}
}

void up_date_z_weight(gsl_matrix *Y, gsl_matrix *YTrans, /*gsl_matrix *WeightedY,*/ gsl_matrix *Mu, gsl_matrix *Precision, gsl_vector *W, gsl_matrix *Z, gsl_matrix *Weight, gsl_vector *SumZ, gsl_vector *SumWZ, double nu, double lambda, double *logLike, int last, int transform)
{
	int i=0,j=0,k=0,ly=(*Y).size1,py=(*Y).size2,K=(*Mu).size1;
	gsl_vector_view row,row1,row2,row3,row4;
	gsl_matrix_view Prow;
	double normConstant=0;
	double like=0,tmpLike=0,logJacobian=0;
	/*	Initialize elements to zero*/
	*logLike=0;
	gsl_vector_set_zero(SumZ);
	gsl_vector_set_zero(SumWZ);

	for(i=0;i<ly;i++)
	{        
		for(j=0;j<py;j++)
		{
			logJacobian+=(lambda-1)*log(fabs(gsl_matrix_get(Y,i,j)));
			gsl_matrix_set(YTrans,i,j,(sgn(gsl_matrix_get(Y,i,j))*pow(fabs(gsl_matrix_get(Y,i,j)),lambda)-1)/lambda);
		}
		row1=gsl_matrix_row(YTrans,i);
		// row3=gsl_matrix_row(WeightedY,i);                 
		normConstant=0;
		like=0;
		for(k=0;k<K;k++)
		{
			row2=gsl_matrix_row(Mu,k);
			row=gsl_matrix_row(Precision,k);
			Prow=gsl_matrix_view_vector(&row.vector,py,py);
			/* Here I assume the cholesky decomposition has been done */			
			tmpLike=gsl_vector_get(W,k)*gsl_ran_mvnt_pdf(&row1.vector,&row2.vector,&Prow.matrix,nu,1,0);            
			/* E[Z|y] */
			gsl_matrix_set(Z,i,k,tmpLike);
			/*  Compute the likelihood */
			like+=tmpLike;
			/* Compute the normlizing constant */
			normConstant+=gsl_matrix_get(Z,i,k);
			/* Compute E[u|y,z] */
			gsl_matrix_set(Weight,i,k,(py+nu)/(gsl_pow_2(gsl_mahalanobis(&Prow.matrix, &row1.vector, &row2.vector, 1))+nu));
		}
		*logLike+=log(like);
		/* Scale the Z's so that they sum to one */
		row=gsl_matrix_row(Z,i);
		gsl_blas_dscal(1.0/normConstant,&row.vector);
		for(k=0;k<K;k++)
		{
			/* Here I assume the cholesky decomposition has been done */
			/* Note that I don't store the weights, I only store WZ */		  
			/* I store the square root for efficiency */
			if(last==0)
			{
				gsl_vector_set(SumZ,k,gsl_vector_get(SumZ,k)+gsl_matrix_get(Z,i,k));
                /* Compute zu */
				gsl_matrix_set(Z,i,k,gsl_matrix_get(Z,i,k)*gsl_matrix_get(Weight,i,k));
				gsl_vector_set(SumWZ,k,gsl_vector_get(SumWZ,k)+gsl_matrix_get(Z,i,k));
                /* Compute (zu)^(1/2) */
				gsl_matrix_set(Weight,i,k,sqrt(gsl_matrix_get(Z,i,k)));
			}
		}
	}
	if(transform == 1) 
	{
		*logLike+=logJacobian;
	}

}

double log_likelihood(gsl_matrix *Y, gsl_matrix *Mu, gsl_matrix *Precision, gsl_vector *W, double nu)
{
	int i=0,k=0;
	int ly=(*Y).size1,py=(*Y).size2, K=(*W).size;
	double log_like=0,like=0;
	gsl_vector_view row,MuRow,YRow;  
	gsl_matrix_view Prow;

	for(i=0;i<ly;i++)
	{
		like=0;
		YRow=gsl_matrix_row(Y,i);
		for(k=0;k<K;k++)
		{	  			
			MuRow=gsl_matrix_row(Mu,k);
			row=gsl_matrix_row(Precision,k);
			Prow=gsl_matrix_view_vector(&row.vector,py,py);
			/* Here I assume the cholesky decomposition has been done */
			like+=gsl_vector_get(W,k)*gsl_ran_mvnt_pdf(&YRow.vector,&MuRow.vector, &Prow.matrix, nu, 1, 0); 

		}
		log_like+=log(like);
	}  
	return(log_like);
}


