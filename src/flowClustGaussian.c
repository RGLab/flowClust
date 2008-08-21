#include "flowClust.h"
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

void flowClustGaussian(double *y, int *ly, int *py, int *K, double *w, double *mu, double *precision, double *lambda, double *z, double *u, int *label, double *uncertainty, double *q_cutoff, double *z_cutoff, int *flagOutliers, int *B, double *tol, int *transform, double *logLike)
{

    gsl_matrix_view Y, Mu, Precision, Z, U;
    gsl_vector_view W;
	gsl_matrix *YTrans=gsl_matrix_alloc(*ly,*py);
	gsl_matrix *ZY=gsl_matrix_calloc(*ly,*K**py);
	gsl_matrix *DiagOne=gsl_matrix_calloc(*py,*py);
	gsl_vector *SumZ=gsl_vector_calloc(*K);
    gsl_vector_view rowMu, rowPrecision, rowYTrans, rowZY, subRowZY, rowZ;
    gsl_matrix_view matrixPrecision, matrixZY;
	int i=0, j=0, k=0;
    double logLikeOld=0;
    double Diff=100.0;  // difference between logLike and logLikeOld
	int iter=0;    // counter for EM iterations
	int BSolve=100;    // Only estimate lambda in the first 100 EM iterations

	// const gsl_rng_type * T;
    // gsl_rng * r;
	const gsl_root_fsolver_type *S;    // the instance of a solver (to be initialized with the desired algorithm, ie, Brent here)
	gsl_root_fsolver *s;    // pointer to a newly allocated instance of a solver
	gsl_function F;
	struct BoxCox_params params;

	double DiffSolve=1e-3;    // tolerance for Brent's algorithm in each M-step
    int iterSolve=0;    // counter for the Brent's algorithm
	int iterSolveMax=50;    // max no. of iterations for the Brent's algorithm in each M-step;
	double lambdaHat=0;    // current estimate of lambda
	double xLow=.1, xUp=1;  // initial search interval of Brent's algorithm
	double x_lo= .1, x_hi=1;    // current bracketing interval for solver
    int status=0;   // status of the solver (an error indicator)    
    
    
	
	/* Create matrix and vector views*/  
	Y=gsl_matrix_view_array(y,*ly,*py);
	W=gsl_vector_view_array(w,*K);
	Mu=gsl_matrix_view_array(mu,*K,*py);
	Precision=gsl_matrix_view_array(precision,*K,*py**py);
	Z=gsl_matrix_view_array(z,*ly,*K);
	U=gsl_matrix_view_array(u,*ly,*K);
	gsl_matrix_set_identity(DiagOne);	
	/* Add a constant to all values to avoid problems with zeros */	
  gsl_matrix_add_constant(&Y.matrix, 0.000001111);
    /* Set xLow to -1 if all observations are positive */
    // if(gsl_matrix_min(&Y.matrix)>0) xLow=-1;

	/*** Initialization ***/

    if (*transform==0) gsl_matrix_memcpy(YTrans, &Y.matrix);
	for(i=0;i<*ly;i++)
	{
        /* Initialize YTrans */
        if (*transform==1)
        {
    		for(j=0;j<*py;j++)
	    	{
		    	gsl_matrix_set(YTrans, i, j, (sgn(gsl_matrix_get(&Y.matrix,i,j)) * pow(fabs(gsl_matrix_get(&Y.matrix,i,j)),*lambda)-1) / (*lambda));
		    }
		}
		/* Initialize Z-matrix (with 1's and 0's) and SumZ according to initial partition */
		if(label[i]>0)
		{
			gsl_matrix_set(&Z.matrix,i,label[i]-1,1.0);
			gsl_vector_set(SumZ,label[i]-1,gsl_vector_get(SumZ,label[i]-1)+1.0);
		}
		/* Initialize ZY-matrix */
		rowZY=gsl_matrix_row(ZY,i);
		rowYTrans=gsl_matrix_row(YTrans,i);        
		for(k=0;k<*K;k++)
		{
			subRowZY=gsl_vector_subvector(&rowZY.vector, k**py, *py);
			gsl_blas_dcopy(&rowYTrans.vector,&subRowZY.vector);
			gsl_blas_dscal(gsl_matrix_get(&Z.matrix,i,k), &subRowZY.vector);
		}
	}

	/* Initialize Mu using BLAS */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans, 0, &Mu.matrix);
	for(k=0;k<*K;k++)
	{
		rowMu=gsl_matrix_row(&Mu.matrix,k);
		gsl_blas_dscal(1./gsl_vector_get(SumZ,k),&rowMu.vector);
		/* Initialize Precision (cluster specific) */
		rowPrecision=gsl_matrix_row(&Precision.matrix,k);
		matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,*py,*py);            
		matrixZY=gsl_matrix_submatrix(ZY, 0, k**py, *ly, *py);      
		up_date_precision(&matrixZY.matrix, &rowMu.vector, &matrixPrecision.matrix, gsl_vector_get(SumZ,k), gsl_vector_get(SumZ,k), DiagOne);            
        /* Initialize Mixing Proportions */
		gsl_vector_set(&W.vector,k,gsl_vector_get(SumZ,k)/(*ly));
	}


	/* This is actually not need but I leave it for now if we want to look at stochastic EM in the future
	gsl_rng_env_setup();
	T=gsl_rng_default;
	r=gsl_rng_alloc(T); */

    /* Initialize the solver */
	S=gsl_root_fsolver_brent;
	s=gsl_root_fsolver_alloc(S);
    /* Initialize the function to solve */
	F.function=&BoxCoxGradient;
	F.params=&params;
	/* Turn off the error handler */
	gsl_set_error_handler_off();

	/* Initialize the log likelihood */
	*logLike=-FLT_MAX;



	/*** EM algorithm ***/
	while((Diff>*tol) & (iter<*B))
	{		
		logLikeOld=*logLike;

		/* E step -- Compute Z, U and also YTrans, SumZ, logLike */
		up_date_z_u_gaussian(&Y.matrix, YTrans, &W.vector, &Mu.matrix, &Precision.matrix, &Z.matrix, &U.matrix, SumZ, *lambda, logLike, *transform, 0);      

		/* M step*/        

		/* Only estimate lambda in the first 100 iterations */
		if(iter>=BSolve) *transform=0;			
			
        /* Solve for lambda */
		if(*transform==1)
		{
    		params.Y=&Y.matrix;
	    	params.W=&W.vector;     // will be updated when calling solver
		    params.Mu=&Mu.matrix;   // will be updated
    		params.Precision=&Precision.matrix;    // will be updated
	    	params.Z=&Z.matrix;
		    params.U=&U.matrix;
    		params.ZUY=ZY;         // will be updated
	    	params.SumZ=SumZ;
		    params.SumZU=SumZ;
    		params.DiagOne=DiagOne;        
	    	params.lambda=*lambda;

    		F.params=&params;

		    /* Initialize the solver s to use the function F and the initial search interval */
			status=gsl_root_fsolver_set(s, &F, xLow, xUp);
            /* Print warning message if error occurs at the initial stage of Brent's algorithm
			if(status!=0)
			{
				warning("Brent's algorithm returned an error, error code %d\n", status);				
			} */

    		if(status==0)
    		{
	    		iterSolve=0;            
		    	do
			    {
				    iterSolve++;
    				/* Perform one iteration of Brent's algorithm */
	    			status=gsl_root_fsolver_iterate(s);
	    			/* current estimate of lambda */
		    		lambdaHat=gsl_root_fsolver_root(s);
    				/* current bracketing interval for solver */
    				x_lo=gsl_root_fsolver_x_lower(s);
	    			x_hi=gsl_root_fsolver_x_upper(s);
		    		status=gsl_root_test_interval(x_lo, x_hi, DiffSolve, 0);
			    }
    			while((status == GSL_CONTINUE) && (iterSolve < iterSolveMax));
                /* Print warning message if error occurs in the Brent's algorithm
		    	if((status!=0) && (iterSolve<iterSolveMax))
			    {
				    warning("Warning in the inner loop, Brent's algorithm returned an error, error code %d\n", status);
    			} */			
                /* Set the new value of the transformation */
		    	*lambda=lambdaHat;        
        	    // Rprintf("%d iterations, lambda = %g\n",iterSolve,lambdaHat);
		    }
		}

        /* Update ZY, Mu, Precision and W */
		if((*transform==0) || ((*transform==1) && (status!=0)))
		{
			/* Update ZY-matrix */
			for(i=0;i<*ly;i++)
			{
				rowYTrans=gsl_matrix_row(YTrans,i);
				rowZY=gsl_matrix_row(ZY,i);                    
				for(k=0;k<*K;k++)
				{                        
					subRowZY=gsl_vector_subvector(&rowZY.vector, k**py, *py);
					gsl_blas_dcopy(&rowYTrans.vector,&subRowZY.vector);                        
					gsl_blas_dscal(gsl_matrix_get(&U.matrix,i,k), &subRowZY.vector);                        
				}
			}

			/* Update Mu using BLAS */
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans, 0, &Mu.matrix);        
			for(k=0;k<*K;k++)
			{
				rowMu=gsl_matrix_row(&Mu.matrix,k);
				gsl_blas_dscal(1./gsl_vector_get(SumZ,k),&rowMu.vector);    
		        /* Update Precision (cluster specific) */
				rowPrecision=gsl_matrix_row(&Precision.matrix,k);
				matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector, *py,*py);            
				matrixZY=gsl_matrix_submatrix(ZY, 0, k**py, *ly, *py);
				up_date_precision(&matrixZY.matrix, &rowMu.vector, &matrixPrecision.matrix, gsl_vector_get(SumZ,k), gsl_vector_get(SumZ,k), DiagOne);                        
                /* Update Mixing Proportions */
				gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ,k)/(*ly));            
			}
		}
		
		Diff=fabs(*logLike-logLikeOld)/fabs(logLikeOld);
		iter++; 
	}
	// Rprintf("The EM required %d iterations\n",iter);
	// Rprintf("The tolerance is %g\n",Diff);

	/* One more E-step to compute the final z's and u's */
	up_date_z_u_gaussian(&Y.matrix, YTrans, &W.vector, &Mu.matrix, &Precision.matrix, &Z.matrix, &U.matrix, SumZ, *lambda, logLike, *transform, 1);      

    /* Output Precision to be the covariance matrix */
	for(k=0;k<*K;k++)
	{
	    gsl_matrix_set_identity(DiagOne);
		rowPrecision=gsl_matrix_row(&Precision.matrix,k);
		matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,*py,*py);
        gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &matrixPrecision.matrix, DiagOne);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, &matrixPrecision.matrix);
    }

	/* Compute uncertainty and flagOutliers */
	for(i=0;i<*ly;i++)
	{
		rowZ=gsl_matrix_row(&Z.matrix,i);
		label[i]=gsl_vector_max_index(&rowZ.vector);
		uncertainty[i]=1-gsl_vector_get(&rowZ.vector,label[i]);
		if((*q_cutoff!=-1 && gsl_matrix_get(&U.matrix,i,label[i])>*q_cutoff) || gsl_matrix_get(&Z.matrix,i,label[i])<*z_cutoff)
			flagOutliers[i]=1;
        /*	else
			flagOutliers[i]=0; */
	}
	
	
	gsl_vector_free(SumZ);
	gsl_matrix_free(DiagOne);
	gsl_matrix_free(ZY);
	gsl_matrix_free(YTrans);
	gsl_root_fsolver_free(s);
}



/* Compute Z, U and also YTrans, SumZ, logLike */
void up_date_z_u_gaussian(gsl_matrix *Y, gsl_matrix *YTrans, gsl_vector *W, gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U, gsl_vector *SumZ, double lambda, double *logLike, int transform, int last)
{
	int i=0,j=0,k=0, ly=(*Y).size1, py=(*Y).size2, K=(*Mu).size1;
	gsl_vector_view rowYTrans, rowMu, rowPrecision, rowZ;
	gsl_matrix_view matrixPrecision;
	double logJacobian=0;
	double normConstant=0;    // normalizing constant for Z
	double tmpLike=0;    // W * density function
	double like=0;    // Sum of (W * density function) wrt one observation
	
	/*	Initialize elements to zero*/
	*logLike=0;
	gsl_vector_set_zero(SumZ);

	for(i=0;i<ly;i++)
	{        
        if(transform == 1) 
        {
    		for(j=0;j<py;j++)
	    	{
		    	logJacobian+=(lambda-1)*log(fabs(gsl_matrix_get(Y,i,j)));
			    gsl_matrix_set(YTrans, i, j, (sgn(gsl_matrix_get(Y,i,j)) * pow(fabs(gsl_matrix_get(Y,i,j)),lambda)-1) / lambda);
		    }
		}
		rowYTrans=gsl_matrix_row(YTrans,i);
		normConstant=0;
		like=0;
		for(k=0;k<K;k++)
		{
			rowMu=gsl_matrix_row(Mu,k);
			rowPrecision=gsl_matrix_row(Precision,k);
			matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,py,py);
			/* Here I assume the cholesky decomposition has been done */
			/* Compute W * density function */
			tmpLike=gsl_vector_get(W,k) * gsl_ran_mvngaussian_pdf(&rowYTrans.vector,&rowMu.vector,&matrixPrecision.matrix,1,0);            
			/* E[Z|y] (un-normalized) */
			gsl_matrix_set(Z,i,k,tmpLike);
			/* Compute the normalizing constant (for Z) */
			normConstant+=gsl_matrix_get(Z,i,k);
			/* Compute U = Mahalanobis distance */
			if(last==1) gsl_matrix_set(U,i,k,gsl_pow_2(gsl_mahalanobis(&matrixPrecision.matrix, &rowYTrans.vector, &rowMu.vector, 1)));
			/*  Compute the likelihood wrt to one observation*/
			like+=tmpLike;
		}
		*logLike+=log(like);
		/* Scale the Z's so that they sum to one */
		rowZ=gsl_matrix_row(Z,i);
		gsl_blas_dscal(1.0/normConstant,&rowZ.vector);

		for(k=0;k<K;k++)
		{
			/* Here I assume the cholesky decomposition has been done */
			/* Note that Z is used to store ZU */
			/* Note that U is used to store (ZU)^(1/2) */						
			/* I store the square root for efficiency */
			if(last==0)    // last=0 means not the final EM iteration
			{
				gsl_vector_set(SumZ, k, gsl_vector_get(SumZ,k) + gsl_matrix_get(Z,i,k));
                /* Compute U=Z^(1/2) */
				gsl_matrix_set(U,i,k,sqrt(gsl_matrix_get(Z,i,k)));
			}
		}
	}
	if(transform == 1) 
	{
		*logLike+=logJacobian;
	}

}



void getEstimatesGaussian(double *y, int *ly, int *py, int *K, double *mu, double *precision, double *z)
{	
	int i=0,j=0,k=0;
	gsl_vector *SumZ=gsl_vector_calloc(*K);
	gsl_matrix *TempMatrix=gsl_matrix_calloc(*py,*py);
	gsl_matrix *DiagOne=gsl_matrix_calloc(*py,*py);
	gsl_matrix *WeightedY=gsl_matrix_calloc(*ly,*py);    // WeightedY = Z*Y
	gsl_matrix_view Y, Mu, Precision, Z, matrixPrecision;
	gsl_vector_view rowMu, rowPrecision, rowWeightedY;
	
	/* Create matrix and vector views*/  
	Y=gsl_matrix_view_array(y,*ly,*py);
	Mu=gsl_matrix_view_array(mu,*K,*py);
	Precision=gsl_matrix_view_array(precision,*K,*py**py);
	Z=gsl_matrix_view_array(z,*ly,*K);
    gsl_matrix_set_identity(DiagOne);

	for(i=0;i<*ly;i++)
	{        
    	for(k=0;k<*K;k++)
	    {
		    gsl_vector_set(SumZ, k, gsl_vector_get(SumZ,k) + gsl_matrix_get(&Z.matrix,i,k));
	    }
    }

	/* Compute Mu using BLAS */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, &Y.matrix, 0.0, &Mu.matrix);        
	for(k=0;k<*K;k++)
	{
		rowMu=gsl_matrix_row(&Mu.matrix,k);
		gsl_blas_dscal(1./gsl_vector_get(SumZ,k),&rowMu.vector);    

        /* Compute WeightedY = Z*Y */
        gsl_matrix_memcpy(WeightedY, &Y.matrix); 
        for(i=0;i<*ly;i++)
        {
            rowWeightedY=gsl_matrix_row(WeightedY,i);
            gsl_blas_dscal(gsl_matrix_get(&Z.matrix,i,k),&rowWeightedY.vector);
        }
        
		/* Compute Precision (cluster specific) */
		rowPrecision=gsl_matrix_row(&Precision.matrix,k);
		matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,*py,*py);            
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0/gsl_vector_get(SumZ,k), WeightedY, &Y.matrix, 0.0, TempMatrix);
    	gsl_blas_dsyr(CblasLower, -1., &rowMu.vector, TempMatrix);    
        gsl_blas_dsymm(CblasLeft, CblasLower, 1.0, TempMatrix, DiagOne, 0.0, &matrixPrecision.matrix);
	}

    gsl_vector_free(SumZ);
    gsl_matrix_free(TempMatrix);
    gsl_matrix_free(DiagOne);
    gsl_matrix_free(WeightedY);
}

