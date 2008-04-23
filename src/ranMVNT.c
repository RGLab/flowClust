#include "flowClust.h"

double gsl_ran_mvngaussian_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, int is_chol, int is_log)
{
	int i=0,n=(*Mu).size;
	double pdf=0;

	gsl_vector *YTilde=gsl_vector_calloc(n);  
	double sum_square=0;
	gsl_matrix *PrecisionTmp=NULL;
	/* Check if the cholesky decomposition has already been computed */
	/* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  		
		PrecisionTmp=gsl_matrix_calloc(n,n);
				/* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
				/* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}

	pdf=-0.5*n*gsl_sf_log(2.*M_PI);

	for(i=0;i<n;i++)
	{
		pdf+=gsl_sf_log(gsl_matrix_get(Precision,i,i));
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

	/* Based on the precision matrix */
	/* Compute (L'Y) */
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
	/* Compute the norm of YTilde=||(L'Y)'LY|| */
	sum_square=gsl_blas_dnrm2(YTilde);
	pdf+=-0.5*gsl_pow_2(sum_square);

	if(is_log==0)
	{
	//pdf=gsl_sf_exp(pdf);
	//gsl_sf_exp_e(pdf, &result);
	//pdf=result.val;
		pdf=exp(pdf);
	}

	if(is_chol==0)
	{
		gsl_matrix_free(PrecisionTmp);
		/* Go back to the original matrix, before chol decomposition */
		gsl_matrix_memcpy(Precision,PrecisionTmp);

	}
	gsl_vector_free(YTilde);


	return(pdf);
}


double gsl_ran_mvnt_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, int is_log)
{
	int i=0,n=(*Mu).size;
	double pdf=0;
	gsl_matrix *PrecisionTmp=NULL;
	gsl_vector *YTilde=gsl_vector_calloc(n);  
	double sum_square=0;

		/* Check if the cholesky decomposition has already been computed */
		/* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(n,n);
				/* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
				/* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}

	for(i=0;i<n;i++)
	{
		pdf+=gsl_sf_log(gsl_matrix_get(Precision,i,i));
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

		/** Based on the precision matrix **/
		/** Compute (L'Y) **/
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
		/** Compute the norm of YTilde=||(L'Y)'LY|| **/
	sum_square=gsl_blas_dnrm2(YTilde);

	pdf+=-0.5*(nu+n)*log(1.0+gsl_pow_2(sum_square)/nu);
	pdf+=gsl_sf_lngamma(0.5*(nu+n))-gsl_sf_lngamma(0.5*nu)-0.5*n*log(nu*M_PI);
	if(is_log==0)
	{
		pdf=exp(pdf);
	}

	if(is_chol==0)
	{
		/* Go back to the original matrix, before chol decomposition   */
		gsl_matrix_memcpy(Precision,PrecisionTmp);
		gsl_matrix_free(PrecisionTmp);
	}

	
	gsl_vector_free(YTilde);

	return(pdf);
}

void gsl_ran_mvngaussian(gsl_vector *Mu, gsl_matrix *Precision, int is_chol, gsl_vector *Y, gsl_rng *r)
{
	int i=0,n=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;

	/* Check if the cholesky decomposition has already been computed */
	/* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(n,n);
		/* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
		/* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}

	for(i=0;i<n;i++)
		gsl_vector_set(Y,i,gsl_ran_gaussian(r,1));

	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, Y);  
	gsl_vector_add(Y,Mu);
	if(is_chol==0)
	{
		/* Go back to the original matrix, before chol decomposition */
		gsl_matrix_memcpy(Precision,PrecisionTmp);
		gsl_matrix_free(PrecisionTmp);
	}

	
}

void gsl_ran_mvnt(gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, gsl_vector *Y, gsl_rng *r)
{
	int i=0,n=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;

		/* Check if the cholesky decomposition has already been computed */
		/* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(n,n);
		/* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
		/* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}
	else
	{
		/** In this case I just set it to zero **/
		PrecisionTmp=gsl_matrix_calloc(1,1);
	}

	for(i=0;i<n;i++)
		gsl_vector_set(Y,i,gsl_ran_gaussian(r,1));
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, Y);  
	gsl_vector_add(Y,Mu);
 	/*	Scale the norm by a Gamma(nu/2,nu/2) to get a multivariate t with nu df */
	gsl_vector_scale(Y,1./gsl_ran_gamma(r,nu/2.,2./nu));
	if(is_chol==0)
	{
		gsl_matrix_free(PrecisionTmp);
		/* Go back to the original matrix, before chol decomposition */
		gsl_matrix_memcpy(Precision,PrecisionTmp);
	}
	

}

double gsl_mahalanobis(gsl_matrix *Precision, gsl_vector *Y, gsl_vector *Mu, int is_chol)
{

	int i=0,n=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;
	gsl_vector *YTilde=gsl_vector_calloc(n);  
	double sum_square=0;

	/* Check if the cholesky decomposition has already been computed */
	/* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(n,n);
				/* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
				/* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}
	else
	{
		/** In this case I just set it to zero **/		
		PrecisionTmp=gsl_matrix_calloc(1,1);
	}

	for(i=0;i<n;i++)
	{
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

	/* Based on the precision matrix */
	/* Compute (L'Y) */
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
	/* Compute the norm of YTilde=||(L'Y)'LY|| */
	sum_square=gsl_blas_dnrm2(YTilde);

	if(is_chol==0)
	{
		gsl_matrix_free(PrecisionTmp);
		/*  Go back to the original matrix, before chol decomposition  */
		gsl_matrix_memcpy(Precision,PrecisionTmp);
	}

	gsl_vector_free(YTilde);
	return(sum_square);
}

