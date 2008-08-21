#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

struct BoxCox_params
{
    gsl_matrix *Y, *Mu, *Precision, *Z, *U, *ZUY, *DiagOne;
    gsl_vector *W, *SumZ, *SumZU;
    double lambda;
};

double BoxCoxGradient(double x, void *params);
int sgn(double x);

void up_date_precision(gsl_matrix *ZUY, gsl_vector *Mu, gsl_matrix *Precision, double SumZ, double SumZU, gsl_matrix *DiagOne);
void up_date_z_u(gsl_matrix *Y, gsl_matrix *YTrans, gsl_vector *W, gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U, gsl_vector *SumZ, gsl_vector *SumZU, double nu, double lambda, double *logLike, int transform, int last);
double log_likelihood(gsl_matrix *Y, gsl_matrix *Mu, gsl_matrix *Precision, gsl_vector *W, double nu);

double gsl_ran_mvngaussian_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, int is_chol, int is_log);
double gsl_ran_mvnt_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, int is_log);
void gsl_ran_mvngaussian(gsl_vector *Mu, gsl_matrix *Precision, int is_chol, gsl_vector *Y, gsl_rng *r);
void gsl_ran_mvnt(gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, gsl_vector *Y, gsl_rng *r);
double gsl_mahalanobis(gsl_matrix *Precision, gsl_vector *Y, gsl_vector *Mu, int is_chol);

void flowClust(double *y, int *ly, int *py, int *K, double *w, double *mu, double *precision, double *lambda, double *nu, double *z, double *u, int *label, double *uncertainty, double *u_cutoff, double *z_cutoff, int *flagOutliers, int *B, double *tol, int *transform, double *logLike);
void getEstimates(double *y, int *ly, int *py, int *K, double *mu, double *precision, double *nu, double *z, double *u);

void flowClustGaussian(double *y, int *ly, int *py, int *K, double *w, double *mu, double *precision, double *lambda, double *z, double *u, int *label, double *uncertainty, double *q_cutoff, double *z_cutoff, int *flagOutliers, int *B, double *tol, int *transform, double *logLike);
void up_date_z_u_gaussian(gsl_matrix *Y, gsl_matrix *YTrans, gsl_vector *W, gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U, gsl_vector *SumZ, double lambda, double *logLike, int transform, int last);
void getEstimatesGaussian(double *y, int *ly, int *py, int *K, double *mu, double *precision, double *z);
