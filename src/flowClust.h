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
    gsl_matrix *Mu, *Precision, *Weight, *Z, *Y, *WeightedY, *DiagOne;
    gsl_vector *SumWZ, *SumZ, *W;
    double lambda;
};

double BoxCoxGradient(double x, void *params);
int sgn(double x);
// double BoxCox(double x, void *params);

// void up_date_mu(gsl_matrix *Y, gsl_vector *Mu, gsl_vector *Weight, double SumWZ);
// void up_date_mu_all(gsl_matrix *Y, gsl_matrix *Mu, gsl_matrix *Z, gsl_matrix *Weight);
void up_date_precision(gsl_matrix *Y, gsl_vector *Mu, gsl_matrix *Precision, /* gsl_vector *Weight, */ double SumZ, double SumWZ, gsl_matrix *DiagOne);
void up_date_z_weight(gsl_matrix *Y, gsl_matrix *YTrans, /*gsl_matrix *WeightedY,*/ gsl_matrix *Mu, gsl_matrix *Precision, gsl_vector *W, gsl_matrix *Z, gsl_matrix *Weight, gsl_vector *SumZ, gsl_vector *SumWZ, double nu, double lambda, double *logLike, int last, int transform);
double log_likelihood(gsl_matrix *Y, gsl_matrix *Mu, gsl_matrix *Precision, gsl_vector *W, double nu);

double gsl_ran_mvngaussian_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, int is_chol, int is_log);
double gsl_ran_mvnt_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, int is_log);
void gsl_ran_mvngaussian(gsl_vector *Mu, gsl_matrix *Precision, int is_chol, gsl_vector *Y, gsl_rng *r);
void gsl_ran_mvnt(gsl_vector *Mu, gsl_matrix *Precision, double nu, int is_chol, gsl_vector *Y, gsl_rng *r);
double gsl_mahalanobis(gsl_matrix *Precision, gsl_vector *Y, gsl_vector *Mu, int is_chol);

void flowClust(double *y, int *ly, int *py, double *mu, double *precision, double *w, double *z, double *u, double *lambda, int *label, double *uncertainty, double *u_cutoff, double *z_cutoff, int *flagOutliers, int *K, double *nu, int *B, double *tol, int *transform, double *logLike);
