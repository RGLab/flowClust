#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <Rinternals.h>
struct BoxCox_params {
	gsl_matrix *Y, *Mu, *Precision, *Z, *U, *logY, *YTrans, *ZUY, *ZUY2,
			*PrecisionUpdate, *MuUpdate, *DiagOne, *Lambda0, *Mu0, *Omega0;
	gsl_vector *W, *SumZ, *SumZU, *SumZlogY;
	int *model, *itersolve;
	double *lambda;
	int K, iter;
	double kappa0;
	double *nu0;
	int	*oorder;
	double* w0;
};

double BoxCoxGradient(double x, void *params);
int sgn(double x);
double BoxCoxGradientPriors(double x, void *params);

struct Nu_params {
	gsl_vector *SumZ, *SumZU, *SumZlogU;
	double *nu;
	int K;
};

double NuGradient(double x, void *params);

struct Nu_ECME_params {
	gsl_matrix *YTrans, *Mu, *Precision, *logY;
	gsl_vector *W;
	double *lambda, *nu;
	int *transform;
};

double NuLikelihood(double x, void *params);

int up_date_precision(gsl_matrix *ZUY, gsl_vector *Mu, gsl_matrix *Precision,
		double SumZ, double SumZU, gsl_matrix *DiagOne, gsl_vector *Mu0,
		gsl_matrix *Lambda0, double kappa0, double nu0);
void up_date_z_u(gsl_matrix *logY, gsl_matrix *YTrans, gsl_vector *W,
		gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U,
		gsl_vector *SumZ, gsl_vector *SumZU, gsl_vector *SumZlogU, double *nu,
		double *lambda, double *logLike, int transform, int nuEstimate,
		int last);
void up_date_z_uS(gsl_matrix *logY, gsl_matrix *YTransS, gsl_vector *W,
		gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U,
		gsl_vector *SumZ, gsl_vector *SumZU, gsl_vector *SumZlogU,
		gsl_vector *SumZlogY, double *nu, double *lambda, double *logLike,
		int nuEstimate, int last);

double gsl_ran_mvngaussian_pdf(gsl_vector *Y, gsl_vector *Mu,
		gsl_matrix *Precision, int is_chol, int is_log);
double gsl_ran_mvnt_pdf(gsl_vector *Y, gsl_vector *Mu, gsl_matrix *Precision,
		double nu, int is_chol, int is_log);
void gsl_ran_mvngaussian(gsl_vector *Mu, gsl_matrix *Precision, int is_chol,
		gsl_vector *Y, gsl_rng *r);
void gsl_ran_mvnt(gsl_vector *Mu, gsl_matrix *Precision, double nu,
		int is_chol, gsl_vector *Y, gsl_rng *r);
double gsl_mahalanobis(gsl_matrix *Precision, gsl_vector *Y, gsl_vector *Mu,
		int is_chol);
void reorderPriors(gsl_matrix *Mu, gsl_matrix *Mu0, gsl_matrix *Lambda0,
		gsl_matrix *Omega0, double* nu0, double* w0, int* oorder);

void flowClust(double *y, int *ly, int *py, int *K, double *w, double *mu,
		double *precision, double *lambda, double *nu, double *z, double *u,
		int *label, double *uncertainty, double *u_cutoff, double *z_cutoff,
		int *flagOutliers, int *B, double *tol, int *transform,
		int *nuEstimate, double *logLike, int *BSolve, int *iterSolveMax,
		double *DiffSolve, double *xLow, double *xUp, double *nuLow,
		double *nuUp, double *mu0, double *kappa0, double *nu0,
		double *lambda0, double *omega0,double *w0, int *model,int *oorder);
void getEstimates(double *y, int *ly, int *py, int *K, double *mu,
		double *precision, double *nu, double *z, double *u);

void flowClustGaussian(double *y, int *ly, int *py, int *K, double *w,
		double *mu, double *precision, double *lambda, double *z, double *u,
		int *label, double *uncertainty, double *q_cutoff, double *z_cutoff,
		int *flagOutliers, int *B, double *tol, int *transform,
		double *logLike, int *BSolve, int *iterSolveMax, double *DiffSolve,
		double *xLow, double *xUp,double *mu0, double *kappa0, double *nu0,
		double *lambda0, double *omega0, int *model);
void up_date_z_u_gaussian(gsl_matrix *logY, gsl_matrix *YTrans, gsl_vector *W,
		gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U,
		gsl_vector *SumZ, double *lambda, double *logLike, int transform,
		int last);
void up_date_z_u_gaussianS(gsl_matrix *logY, gsl_matrix *YTransS,
		gsl_vector *W, gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z,
		gsl_matrix *U, gsl_vector *SumZ, gsl_vector *SumZlogY, double *lambda,
		double *logLike, int last);
void getEstimatesGaussian(double *y, int *ly, int *py, int *K, double *mu,
		double *precision, double *z);
int ECMUpdateMUg(gsl_matrix *Precision, gsl_vector *Mu, gsl_vector *ZY,
		gsl_matrix *Omega0, double SumZ, gsl_matrix* DiagOne);
int ECMUpdateSigmaG2(gsl_vector* U, gsl_matrix* YTrans, gsl_matrix* MuGMatrix,
		gsl_matrix* Precision, double SumZ, double nu0);
int up_date_precision_gaussian(gsl_matrix *ZUY, gsl_vector *Mu,
		gsl_matrix *Precision, double SumZ, double SumZU, gsl_matrix *DiagOne);
