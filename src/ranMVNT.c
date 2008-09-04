#include "flowClust.h"

/* Compute the density value of the multivariate Gaussian distribution */
double gsl_ran_mvngaussian_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, int is_chol, int is_log)
{
	int i=0, py=(*Mu).size;
	double pdf=0;
	gsl_matrix *PrecisionTmp=NULL;
	gsl_vector *YTilde=gsl_vector_calloc(py);    // Y-Mu
	double sum_square=0;   // store sqrt of Mahalanobis distance

    /* Check if the cholesky decomposition has already been computed */
    /* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(py,py);
        /* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
        /* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}

	pdf=-0.5*py*gsl_sf_log(2.*M_PI);

	for(i=0;i<py;i++)
	{
        /* log of determinant of sqrt of Precision matrix */
		pdf+=gsl_sf_log(gsl_matrix_get(Precision,i,i));
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

    /** Based on the precision matrix **/
    /** Compute (L'Y) **/
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
    /** Compute the norm of YTilde=||(L'Y)'LY|| **/
	sum_square=gsl_blas_dnrm2(YTilde);
	pdf+=-0.5*gsl_pow_2(sum_square);

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


/* Compute the density value of the multivariate t distribution */
double gsl_ran_mvnt_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, int is_log)
{
	int i=0, py=(*Mu).size;
	double pdf=0;
	gsl_matrix *PrecisionTmp=NULL;
	gsl_vector *YTilde=gsl_vector_calloc(py);    // Y-Mu
	double sum_square=0;   // store sqrt of Mahalanobis distance

    /* Check if the cholesky decomposition has already been computed */
    /* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(py,py);
        /* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
        /* Perform the cholesky decomposition*/
		gsl_linalg_cholesky_decomp(Precision);
	}

	for(i=0;i<py;i++)
	{
        /* log of determinant of sqrt of Precision matrix */
		pdf+=gsl_sf_log(gsl_matrix_get(Precision,i,i));
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

    /** Based on the precision matrix **/
    /** Compute (L'Y) **/
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
    /** Compute the norm of YTilde=||(L'Y)'LY|| **/
	sum_square=gsl_blas_dnrm2(YTilde);
	
	pdf+=-0.5*(nu+py)*log(1.0+gsl_pow_2(sum_square)/nu);
    pdf+=gsl_sf_lngamma(0.5*(nu+py))-gsl_sf_lngamma(0.5*nu)-0.5*py*log(nu*M_PI);
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


/* Generate a multivariate Gaussian vector */
void gsl_ran_mvngaussian(gsl_vector *Mu, gsl_matrix *Precision, int is_chol, gsl_vector *Y, gsl_rng *r)
{
	int i=0, py=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;

    /* Check if the cholesky decomposition has already been computed */
    /* This is a bit inefficient and you should do the decomposition before */
	if(is_chol==0)
	{	  
		PrecisionTmp=gsl_matrix_calloc(py,py);
        /* Store the original Precision matrix */
		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
        /* Perform the cholesky decomposition */
		gsl_linalg_cholesky_decomp(Precision);
	}

	for(i=0;i<py;i++)
		gsl_vector_set(Y,i,gsl_ran_gaussian(r,1));
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, Y);  
	gsl_vector_add(Y,Mu);
		
	if(is_chol==0)
  	{
        /* Go back to the original matrix, before chol decomposition   */
  		gsl_matrix_memcpy(Precision,PrecisionTmp);
  		gsl_matrix_free(PrecisionTmp);
  	}
}


/* Generate a multivariate t vector */
void gsl_ran_mvnt(gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, gsl_vector *Y, gsl_rng *r)
{
	int i=0, py=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;

    /* Check if the cholesky decomposition has already been computed */
    /* This is a bit inefficient and you should do the decomposition before */
  	if(is_chol==0)
  	{	  
		PrecisionTmp=gsl_matrix_calloc(py,py);
        /* Store the original Precision matrix */
  		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
        /* Perform the cholesky decomposition */
  		gsl_linalg_cholesky_decomp(Precision);
  	}
  	
	for(i=0;i<py;i++)
		gsl_vector_set(Y,i,gsl_ran_gaussian(r,1));
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, Y);  
    /* Scale the Gaussian vector by a Gamma(nu/2,nu/2) to get a multivariate t vector with nu df */
	gsl_vector_scale(Y,1./sqrt(gsl_ran_gamma(r,nu/2.,2./nu)));
	gsl_vector_add(Y,Mu);

	if(is_chol==0)
  	{
        /* Go back to the original matrix, before chol decomposition */
  		gsl_matrix_memcpy(Precision,PrecisionTmp);
  		gsl_matrix_free(PrecisionTmp);
  	}
}


/* Compute (square root of) Mahalanobis distance */
double gsl_mahalanobis(gsl_matrix *Precision, gsl_vector *Y, gsl_vector *Mu, int is_chol)
{
	int i=0, py=(*Mu).size;
	gsl_matrix *PrecisionTmp=NULL;
	gsl_vector *YTilde=gsl_vector_calloc(py);    // Y-Mu
	double sum_square=0;   // store sqrt of Mahalanobis distance

    /* Check if the cholesky decomposition has already been computed */
    /* This is a bit inefficient and you should do the decomposition before */
  	if(is_chol==0)
  	{	  
		PrecisionTmp=gsl_matrix_calloc(py,py);
        /* Store the original Precision matrix */
  		gsl_matrix_memcpy(PrecisionTmp,Precision);  	
        /* Perform the cholesky decomposition*/
  		gsl_linalg_cholesky_decomp(Precision);
  	}

	for(i=0;i<py;i++)
	{
		gsl_vector_set(YTilde,i,gsl_vector_get(Y,i)-gsl_vector_get(Mu,i));
	}

    /** Based on the precision matrix **/
    /** Compute (L'Y) **/
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, Precision, YTilde);  
    /** Compute the norm of YTilde=||(L'Y)'LY|| **/
	sum_square=gsl_blas_dnrm2(YTilde);

	if(is_chol==0)
  	{
        /*  Go back to the original matrix, before chol decomposition */
  		gsl_matrix_memcpy(Precision,PrecisionTmp);
  		gsl_matrix_free(PrecisionTmp);
  	}

	gsl_vector_free(YTilde);

	return(sum_square);
}

