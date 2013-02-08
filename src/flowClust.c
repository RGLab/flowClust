#include "flowClust.h"
#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>

static const R_CMethodDef CEntries[] = { { "flowClust", (DL_FUNC) & flowClust,
		20 }, { NULL, NULL, 0 } };

// intended to comment out R_init_flowClust. We do not encourage users to call the C routines in R.  We have already defined the corresponding R functions for users to use.
/* void R_init_flowClust(DllInfo *dll)
 {
 R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
 // R_useDynamicSymbols(dll, FALSE);
 } */

void flowClust(double *y, int *ly, int *py, int *K, double *w, double *mu,
		double *precision, double *lambda, double *nu, double *z, double *u,
		int *label, double *uncertainty, double *u_cutoff, double *z_cutoff,
		int *flagOutliers, int *B, double *tol, int *transform,
		int *nuEstimate, double *logLike, int *BSolve, int *iterSolveMax,
		double *DiffSolve, double *xLow, double *xUp, double *nuLow,
		double *nuUp, double *mu0, double *kappa0, double *nu0,
		double *lambda0, double *omega0, double* w0,int *model,  int *oorder) {
			
	gsl_matrix_view Y, Mu, Precision, Mu0, Lambda0, Omega0;
	gsl_matrix_view Z; // zu at intermediate steps; z at output
	gsl_matrix_view U; // sqrt(zu) at intermediate steps; u at output
	gsl_vector_view W;
	gsl_matrix *logY = NULL; // log(abs(Y))
	gsl_matrix *tmpOmega = gsl_matrix_calloc(*K,*py * *py); //working variable for calculations on mean
	gsl_matrix *YTrans = gsl_matrix_calloc(*ly, *py);   
	gsl_matrix *MuGMatrix = gsl_matrix_calloc(*ly, *py);
	gsl_matrix *YTransS = NULL;
	gsl_matrix *ZY = gsl_matrix_calloc(*K, *py); // ZU*YTrans
	gsl_matrix *ZUY = gsl_matrix_calloc(*ly, *K * *py); // sqrt(zu) * YTrans
	gsl_matrix *DiagOne = gsl_matrix_calloc(*py, *py);
	gsl_vector *SumZ = gsl_vector_calloc(*K); // sum(z)
	gsl_vector *SumZU = gsl_vector_calloc(*K); // sum(zu)
	gsl_vector *SumZlogU = NULL; // sum(z*log(u))
	gsl_vector *SumZlogY = NULL; // sum(z*log(abs(y)))
	gsl_vector_view rowMu, rowMu0, rowPrecision, rowLambda0, rowYTrans,
	rowYTransS, rowZUY, subRowZUY, rowZ, rowZUY2, subRowZUY2;
	gsl_matrix_view matrixPrecision, matrixLambda0, matrixYTransS, matrixZUY,
	matrixZUY2;

	int i = 0, j = 0, k = 0;
	double logLikeOld = 0;
	double Diff = 100.0; // difference between logLike and logLikeOld
	int iter = 0; // counter for EM iterations
	// int BSolve=100;    // Only estimate lambda in the first 100 EM iterations

	// const gsl_rng_type * T;
	// gsl_rng * r;
	const gsl_root_fsolver_type *S; // the instance of a solver (to be initialized with the desired algorithm, ie, Brent here)
	gsl_root_fsolver *s, *s2; // pointer to a newly allocated instance of a solver
	const gsl_min_fminimizer_type *S3; // the instance of a minimizer (to be initialized with the desired algorith, ie, Brent here)
	gsl_min_fminimizer *s3; // pointer to a newly allocated instance of a minimizer
	gsl_function F, F2, F3;
	struct BoxCox_params params;
	struct Nu_params params2;
	struct Nu_ECME_params params3;

	// double DiffSolve=1e-3;    // tolerance for Brent's algorithm in each M-step
	int iterSolve = 0; // counter for the Brent's algorithm
	// int iterSolveMax=50;    // max no. of iterations for the Brent's algorithm in each M-step;
	double lambdaHat = 0, nuHat = 0; // current estimate of lambda and nu
	// double xLow=.1, xUp=1, nuLow=2, nuUp=30;  // initial search interval of Brent's algorithm
	double x_lo = 0, x_hi = 0, nu_lo = 0, nu_hi = 0; // current bracketing interval for solver
	int status = 0, status2 = 0; // status of the solver (an error indicator)

	/* Turn off the error handler */
	gsl_set_error_handler_off();

	/* Create matrix and vector views*/
	Y = gsl_matrix_view_array(y, *ly, *py);
	W = gsl_vector_view_array(w, *K);
	Mu = gsl_matrix_view_array(mu, *K, *py);
	Precision = gsl_matrix_view_array(precision, *K, *py * *py);
	Mu0 = gsl_matrix_view_array(mu0, *K, *py);
	Lambda0 = gsl_matrix_view_array(lambda0, *K, *py * *py);
	Omega0 = gsl_matrix_view_array(omega0, *K,*py * *py);
	Z = gsl_matrix_view_array(z, *ly, *K);
	U = gsl_matrix_view_array(u, *ly, *K);
	gsl_matrix_set_identity(DiagOne);
	
	/*
		For Dirichlet prior on component weights w.
		w0 are parameters of the Dirichlet.
		Calculate sum w0
	*/
	double w0sum=0;
	for(int ii=0;ii<*K;ii++){
		w0sum=w0sum+w0[ii];
	}
	
	
	if (*transform >= 1 || (*transform == 0 && *lambda != 1))
		logY = gsl_matrix_calloc(*ly, *py);

	if (*transform > 1) {
		YTransS = gsl_matrix_calloc(*ly, *K * *py);
		SumZlogY = gsl_vector_calloc(*K);
	}
	if (*nuEstimate > 0)
		SumZlogU = gsl_vector_calloc(*K);
	/* Add a constant to all values to avoid problems with zeros */
	gsl_matrix_add_constant(&Y.matrix, 0.000001111);
	/* Set xLow to -1 if all observations are positive */
	// if(gsl_matrix_min(&Y.matrix)>0) xLow=-1;

	/*** Initialization ***/

	if (*transform == 0 && *lambda == 1)
		gsl_matrix_memcpy(YTrans, &Y.matrix);
	for (i = 0; i < *ly; i++) {
		/* Initialize YTrans */
		if (*transform >= 1 || (*transform == 0 && *lambda != 1)) {
			for (j = 0; j < *py; j++) {
				gsl_matrix_set(logY, i, j, gsl_sf_log(fabs(gsl_matrix_get(
						&Y.matrix, i, j))));
				gsl_matrix_set(YTrans, i, j, (sgn(gsl_matrix_get(&Y.matrix, i,
						j)) * pow(fabs(gsl_matrix_get(&Y.matrix, i, j)),
								*lambda) - 1) / (*lambda));
			}
		}
		/* Initialize Z-matrix (with 1's and 0's) and SumZ according to initial partition */
		if (label[i] > 0) {
			gsl_matrix_set(&Z.matrix, i, label[i] - 1, 1.0);
			gsl_vector_set(SumZ, label[i] - 1, gsl_vector_get(SumZ, label[i]
			                                                              - 1) + 1.0);
		}
		/* Initialize ZUY-matrix and ZUY2-matrix*/
		rowZUY = gsl_matrix_row(ZUY, i);
		rowYTrans = gsl_matrix_row(YTrans, i);
		for (k = 0; k < *K; k++) {
			subRowZUY = gsl_vector_subvector(&rowZUY.vector, k * *py, *py);
			gsl_blas_dcopy(&rowYTrans.vector, &subRowZUY.vector);
			gsl_blas_dscal(gsl_matrix_get(&Z.matrix, i, k), &subRowZUY.vector);
		}
	}
	gsl_blas_dcopy(SumZ, SumZU);
	if (*transform > 1) {
		for (k = 0; k < *K; k++) {
			matrixYTransS = gsl_matrix_submatrix(YTransS, 0, k * *py, *ly, *py);
			gsl_matrix_memcpy(&matrixYTransS.matrix, YTrans);
		}
	}

	/* initialization code for ECM with non-conjugate priors */
	if (*model == 2) {
		// Initialize Mu to the prior means
		gsl_matrix_memcpy(&Mu.matrix, &Mu0.matrix);
		/*
			Expectation conditional mode algorithm.
			iter % 2, estimate Mug, SigmaG on alternate iterations due to dependance.
		*/
		if ((iter % 2) == 0 || (iter == 0)) {

			gsl_matrix_memcpy(tmpOmega, &Omega0.matrix);

			/* ZY = Z*Y */// dim(k*p)
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans,
					0.0, ZY);

			/* Cluster Specific Updates of Mu*/

			for (k = 0; k < *K; k++) {
				/*Assume Omega0 is already inverted*/
				/*  Mu = Omega0^{-1} * Mu0*/
				gsl_vector_view rowOmega0 = gsl_matrix_row(tmpOmega,k);
				gsl_matrix_view matrixOmega0 = gsl_matrix_view_vector(&rowOmega0.vector,*py,*py);
				rowMu = gsl_matrix_row(&Mu.matrix, k);
				rowMu0=gsl_matrix_row(&Mu0.matrix,k);

				gsl_blas_dgemv(CblasNoTrans, 1.0,
						&matrixOmega0.matrix, &rowMu0.vector, 0.0, &rowMu.vector);

				rowPrecision = gsl_matrix_row(&Precision.matrix, k);
				matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector,
						*py, *py);

				gsl_vector_view rowZY = gsl_matrix_row(ZY, k);
				ECMUpdateMUg(&matrixPrecision.matrix, &rowMu.vector,
						&rowZY.vector, &matrixOmega0.matrix, gsl_vector_get(SumZU, k),
						DiagOne);
				gsl_vector_set(&W.vector, k, (gsl_vector_get(SumZ, k)+w0[k]) / (w0sum+*ly));
			}
		}
		if ((iter % 2 != 0) || (iter == 0)) {
			for (k = 0; k < *K; k++) {
				rowPrecision = gsl_matrix_row(&Precision.matrix, k);
				matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector,
						*py, *py);
				rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
				matrixLambda0 = gsl_matrix_view_vector(&rowLambda0.vector, *py,
						*py);
				matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly, *py);
				rowMu = gsl_matrix_row(&Mu.matrix, k);
				if (gsl_matrix_isnull(MuGMatrix)) {
					for (i = 0; i < *ly; i++) {
						gsl_matrix_set_row(MuGMatrix, i, &rowMu.vector);
					}
				}
				gsl_matrix_memcpy(&matrixPrecision.matrix,
						&matrixLambda0.matrix);
				gsl_vector_view zu = gsl_matrix_column(&Z.matrix,k);
				status2 = ECMUpdateSigmaG2(&zu.vector,YTrans,MuGMatrix,
						&matrixPrecision.matrix, gsl_vector_get(SumZ, k),nu0[k]);
				// Update weights, include priors.
				gsl_vector_set(&W.vector, k, (gsl_vector_get(SumZ, k)+w0[k]) / (w0sum+*ly));
				gsl_matrix_set_zero(MuGMatrix);
			}
		}
	} else {
		/* initialization for standard flowClust */

		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix, YTrans,
				*kappa0, &Mu.matrix);

		for (k = 0; k < *K; k++) {
			rowMu = gsl_matrix_row(&Mu.matrix, k);
			rowMu0 = gsl_matrix_row(&Mu0.matrix, k);
			gsl_blas_dscal(1. / (gsl_vector_get(SumZU, k) + *kappa0),
					&rowMu.vector);
			/* Initialize Precision (cluster specific) */
			rowPrecision = gsl_matrix_row(&Precision.matrix, k);
			rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
			matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector, *py,
					*py);
			matrixLambda0
			= gsl_matrix_view_vector(&rowLambda0.vector, *py, *py);
			matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly, *py);
			status2 = up_date_precision(&matrixZUY.matrix, &rowMu.vector,
					&matrixPrecision.matrix, gsl_vector_get(SumZ, k),
					gsl_vector_get(SumZU, k), DiagOne, &rowMu0.vector,
					&matrixLambda0.matrix, *kappa0, nu0[k]);
			/* Initialize Mixing Proportions */
			gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ, k) / (*ly));
		}
	}
	/* Initialize the solver */
	S = gsl_root_fsolver_brent;
	if (*transform >= 1 && *model == 1) {
		s = gsl_root_fsolver_alloc(S);
		/* Initialize the function to solve */
		F.function = &BoxCoxGradient;
		params.Y = &Y.matrix;
		params.logY = logY;
	}
	if (*transform == 1 && *model == 2) {
		s = gsl_root_fsolver_alloc(S);
		/* Initialize the function to solve */
		F.function = &BoxCoxGradientPriors;
		params.Y = &Y.matrix;
		params.logY = logY;
		params.MuUpdate= gsl_matrix_calloc(*K,*py);
		params.PrecisionUpdate = gsl_matrix_calloc(*K, *py * *py);
	}
	if (*nuEstimate >= 1) {
		s2 = gsl_root_fsolver_alloc(S);
		/* Initialize the function to solve */
		F2.function = &NuGradient;
	}
	S3 = gsl_min_fminimizer_brent;
	if (*nuEstimate == -1) {
		s3 = gsl_min_fminimizer_alloc(S3);
		/* Initialize the function to solve */
		F3.function = &NuLikelihood;
		params3.logY = logY;
		params3.transform = transform;
		if (*transform == 0) {
			params3.YTrans = YTrans;
			params3.lambda = lambda;
		}
	}

	/* Initialize the log likelihood */
	*logLike = -FLT_MAX;

	/**  algorithm **/
	while ((Diff > *tol) & (iter < *B) & (status2 == 0)) {
		/* This allows us to interupt the code with Ctrl-C */
		R_CheckUserInterrupt();
		logLikeOld = *logLike;

		/* E step -- Compute Z, U, (SumZlogU, SumZlogY) and also SumZ, SumZU, logLike */
		if (*transform <= 1) {
			up_date_z_u(logY, YTrans, &W.vector, &Mu.matrix, &Precision.matrix,
					&Z.matrix, &U.matrix, SumZ, SumZU, SumZlogU, nu, lambda,
					logLike, *transform, *nuEstimate, 0);
		} else {
			up_date_z_uS(logY, YTransS, &W.vector, &Mu.matrix,
					&Precision.matrix, &Z.matrix, &U.matrix, SumZ, SumZU,
					SumZlogU, SumZlogY, nu, lambda, logLike, *nuEstimate, 0);
		}

		/* M step*/

		/* Solve for lambda */
		if ((*transform >= 1) && (iter < *BSolve)) // Only estimate lambda in the first 100(BSolve) iterations
		{
			params.W = &W.vector; // will be updated when calling solver
			params.Mu = &Mu.matrix; // will be updated using MuUpdate
			params.Precision = &Precision.matrix; // will be updated using PrecisionUpdate
			params.Mu0 = &Mu0.matrix; // prior parameters
			params.Lambda0 = &Lambda0.matrix; // prior parameters
			params.kappa0 = *kappa0; // prior parameters
			params.nu0 = nu0; // prior parameters
			params.oorder = oorder;
			params.w0 = w0;

			params.Z = &Z.matrix;
			params.U = &U.matrix;
			params.ZUY = ZUY; // will be updated
			params.SumZ = SumZ;
			params.SumZU = SumZU;
			params.SumZlogY = SumZlogY;
			params.DiagOne = DiagOne;
			params.model = model;
			params.iter = iter;
			params.Omega0 = &Omega0.matrix;
			params.itersolve = &iterSolve;

			if (*transform == 1) {
				params.YTrans = YTrans; // will be updated
				params.K = -1;
				F.params = &params;

				/* Initialize the solver s to use the function F and the initial search interval */
				status = gsl_root_fsolver_set(s, &F, *xLow, *xUp);
				//Rprintf("status=%d  %s\n",status,gsl_strerror(status));

				if (status == 0) {
					iterSolve = 0;
					do {
						iterSolve++;
						/* Perform one iteration of Brent's algorithm */
						status = gsl_root_fsolver_iterate(s);
						/* current estimate of lambda */
						lambdaHat = gsl_root_fsolver_root(s);

						/* current bracketing interval for solver */
						x_lo = gsl_root_fsolver_x_lower(s);
						x_hi = gsl_root_fsolver_x_upper(s);
						status = gsl_root_test_interval(x_lo, x_hi, *DiffSolve,
								0);
					} while ((status == GSL_CONTINUE) && (iterSolve
							< *iterSolveMax));
//					/* Set the new value of the transformation */
					*lambda = lambdaHat;

					if(*model==2){
						gsl_matrix_memcpy(&Mu.matrix,params.MuUpdate);
						gsl_matrix_memcpy(&Precision.matrix,params.PrecisionUpdate);
						reorderPriors(&Mu.matrix,&Mu0.matrix,&Lambda0.matrix,&Omega0.matrix,nu0,w0,oorder);
					}
				}

			} else { //else loop for transform == 2
				for (k = 0; k < *K; k++) {
					matrixYTransS = gsl_matrix_submatrix(YTransS, 0, k * *py,
							*ly, *py);
					params.YTrans = &matrixYTransS.matrix; // will be updated
					params.K = k;
					F.params = &params;

					/* Initialize the solver s to use the function F and the initial search interval */
					status = gsl_root_fsolver_set(s, &F, *xLow, *xUp);
					if (status == 0) {
						iterSolve = 0;
						do {
							iterSolve++;
							/* Perform one iteration of Brent's algorithm */
							status = gsl_root_fsolver_iterate(s);
							/* current estimate of lambda */
							lambdaHat = gsl_root_fsolver_root(s);
							/* current bracketing interval for solver */
							x_lo = gsl_root_fsolver_x_lower(s);
							x_hi = gsl_root_fsolver_x_upper(s);
							status = gsl_root_test_interval(x_lo, x_hi,
									*DiffSolve, 0);
						} while ((status == GSL_CONTINUE) && (iterSolve
								< *iterSolveMax));
						/* Set the new value of the transformation */
						lambda[k] = lambdaHat;
					}

					/* Update YTransS, ZUY, Mu, Precision and W (wrt kth cluster) when Brent incurs an error */
					if (status != 0) {
						/* Update YTransS & ZUY-matrix */
						matrixYTransS = gsl_matrix_submatrix(YTransS, 0, k
								* *py, *ly, *py);
						matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly,
								*py);
						for (i = 0; i < *ly; i++) {
							/* Update YTransS */
							rowYTransS = gsl_matrix_row(&matrixYTransS.matrix,
									i);
							for (j = 0; j < *py; j++) {
								gsl_vector_set(&rowYTransS.vector, j, (sgn(
										gsl_matrix_get(&Y.matrix, i, j)) * pow(
												fabs(gsl_matrix_get(&Y.matrix, i, j)),
												lambda[k]) - 1) / (lambda[k]));
							}

							/* Update ZUY-matrix */
							rowZUY = gsl_matrix_row(&matrixZUY.matrix, i);
							gsl_blas_dcopy(&rowYTransS.vector, &rowZUY.vector);
							gsl_blas_dscal(gsl_matrix_get(&U.matrix, i, k),
									&rowZUY.vector);
						}

						/* Update Mu using BLAS */
						rowMu = gsl_matrix_row(&Mu.matrix, k);
						rowZ = gsl_matrix_column(&Z.matrix, k);

						/* Initialize Mu to Mu0 */
						rowMu0 = gsl_matrix_row(&Mu0.matrix, k);
						gsl_vector_memcpy(&rowMu.vector, &rowMu0.vector);
						gsl_blas_dgemv(CblasTrans, 1.0, &matrixYTransS.matrix,
								&rowZ.vector, *kappa0, &rowMu.vector);
						gsl_blas_dscal(1.
								/ (gsl_vector_get(SumZU, k) + *kappa0),
								&rowMu.vector);

						/* Update Precision (cluster specific) */
						rowPrecision = gsl_matrix_row(&Precision.matrix, k);
						matrixPrecision = gsl_matrix_view_vector(
								&rowPrecision.vector, *py, *py);
						rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
						matrixLambda0 = gsl_matrix_view_vector(
								&rowLambda0.vector, *py, *py);
						status2 = up_date_precision(&matrixZUY.matrix,
								&rowMu.vector, &matrixPrecision.matrix,
								gsl_vector_get(SumZ, k), gsl_vector_get(SumZU,
										k), DiagOne, &rowMu0.vector,
										&matrixLambda0.matrix, *kappa0, nu0[k]);

						/* Update Mixing Proportions */
						gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ, k)
								/ (*ly));
					}
				}
			}
		}

		/* Update YTrans, ZUY, Mu, Precision and W */
		if ((*transform == 0) || ((*transform == 1) && (iter >= *BSolve))
				|| ((*transform == 1) && (status != 0))) {
			/* Update YTrans */
			if ((*transform == 1) && (status != 0)) {
				for (i = 0; i < *ly; i++) {
					for (j = 0; j < *py; j++) {
						gsl_matrix_set(YTrans, i, j, (sgn(gsl_matrix_get(
								&Y.matrix, i, j)) * pow(fabs(gsl_matrix_get(
										&Y.matrix, i, j)), *lambda) - 1) / (*lambda));
					}
				}
			}

			/* Update ZUY-matrix */
			for (i = 0; i < *ly; i++) {
				rowYTrans = gsl_matrix_row(YTrans, i);
				rowZUY = gsl_matrix_row(ZUY, i);
				for (k = 0; k < *K; k++) {
					subRowZUY = gsl_vector_subvector(&rowZUY.vector, k * *py,
							*py);
					gsl_blas_dcopy(&rowYTrans.vector, &subRowZUY.vector);
					gsl_blas_dscal(gsl_matrix_get(&U.matrix, i, k),
							&subRowZUY.vector);
				}
			}

			/* Classic flowClust or with a vague prior */
			if (*model != 2) {
				/* Update Mu using BLAS */
				/* Initialize Mu to the prior */
				gsl_matrix_memcpy(&Mu.matrix, &Mu0.matrix);
				gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix,
						YTrans, *kappa0, &Mu.matrix);
				for (k = 0; k < *K; k++) {
					rowMu = gsl_matrix_row(&Mu.matrix, k);
					rowMu0 = gsl_matrix_row(&Mu0.matrix, k);
					gsl_blas_dscal(1. / (gsl_vector_get(SumZU, k) + *kappa0),
							&rowMu.vector);

					/* Update Precision (cluster specific) */
					rowPrecision = gsl_matrix_row(&Precision.matrix, k);
					matrixPrecision = gsl_matrix_view_vector(
							&rowPrecision.vector, *py, *py);
					rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
					matrixLambda0 = gsl_matrix_view_vector(&rowLambda0.vector,
							*py, *py);
					matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly, *py);
					status2 = up_date_precision(&matrixZUY.matrix,
							&rowMu.vector, &matrixPrecision.matrix,
							gsl_vector_get(SumZ, k), gsl_vector_get(SumZU, k),
							DiagOne, &rowMu0.vector, &matrixLambda0.matrix,
							*kappa0, nu0[k]);
					/* Update Mixing Proportions */
					/* no prior */
					gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ, k)
							/ (*ly));
				}
			}
			/* Non-conjugate prior, ECM algorithm. */
			else if (*model == 2) {
				//Update Mu and Sigma using ECM.//
				// Alternate iterations //
				if ((iter % 2) == 0) {
					gsl_matrix_memcpy(tmpOmega, &Omega0.matrix);

					/*Assume Omega0 is already inverted in flowClust.R*/
					/* ZY = Z*Y */// dim(k*p)
					gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Z.matrix,
							YTrans, 0.0, ZY);

					/* Cluster Specific Updates */
					for (k = 0; k < *K; k++) {
						gsl_vector_view rowOmega0 = gsl_matrix_row(tmpOmega,k);
						gsl_matrix_view matrixOmega0 = gsl_matrix_view_vector(&rowOmega0.vector,*py,*py);
						rowMu = gsl_matrix_row(&Mu.matrix, k);
						rowMu0=gsl_matrix_row(&Mu0.matrix,k);

						gsl_blas_dgemv(CblasNoTrans, 1.0,
								&matrixOmega0.matrix, &rowMu0.vector, 0.0, &rowMu.vector);

						rowPrecision = gsl_matrix_row(&Precision.matrix, k);
						matrixPrecision = gsl_matrix_view_vector(
								&rowPrecision.vector, *py, *py);

						gsl_vector_view rowZY = gsl_matrix_row(ZY, k);
						
						/*	Update MuG */
						ECMUpdateMUg(&matrixPrecision.matrix, &rowMu.vector,
								&rowZY.vector, &matrixOmega0.matrix,
								gsl_vector_get(SumZU, k), DiagOne);
						/* update mixing proportions */
						// gsl_vector_set(&W.vector, k, (gsl_vector_get(SumZ, k))
						// 		/ ((*ly)));
						gsl_vector_set(&W.vector, k, (gsl_vector_get(SumZ, k)+w0[k])
								/ (w0sum+(*ly)));
					}
					reorderPriors(&Mu.matrix,&Mu0.matrix,&Lambda0.matrix,&Omega0.matrix,nu0,w0,oorder);
				}  if (iter % 2 != 0) {
					for (k = 0; k < *K; k++) {
						rowPrecision = gsl_matrix_row(&Precision.matrix, k);
						matrixPrecision = gsl_matrix_view_vector(
								&rowPrecision.vector, *py, *py);
						rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
						matrixLambda0 = gsl_matrix_view_vector(
								&rowLambda0.vector, *py, *py);
						matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly,
								*py);
						rowMu = gsl_matrix_row(&Mu.matrix, k);
						if (gsl_matrix_isnull(MuGMatrix)) {
							for (i = 0; i < *ly; i++) {
								gsl_matrix_set_row(MuGMatrix, i, &rowMu.vector);
							}
						}
						gsl_matrix_memcpy(&matrixPrecision.matrix,
								&matrixLambda0.matrix);
						gsl_vector_view vu = gsl_matrix_column(&U.matrix,k);
						status2 = ECMUpdateSigmaG2(&vu.vector,YTrans,MuGMatrix,
								&matrixPrecision.matrix, gsl_vector_get(SumZ, k),nu0[k]);

						/* update mixing proportions */
						// gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ, k)
						// 		/ (*ly));
						gsl_vector_set(&W.vector, k, (gsl_vector_get(SumZ, k)+w0[k])
									/ (w0sum+(*ly)));
								
						gsl_matrix_set_zero(MuGMatrix);
					}
				}
			}
		} else if (*transform > 1 && iter >= *BSolve) {
			/* No support for priors using transform > 1 */
			for (k = 0; k < *K; k++) {
				/* Update ZUY-matrix */
				matrixYTransS = gsl_matrix_submatrix(YTransS, 0, k * *py, *ly,
						*py);
				matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * *py, *ly, *py);
				for (i = 0; i < *ly; i++) {
					rowYTransS = gsl_matrix_row(&matrixYTransS.matrix, i);
					rowZUY = gsl_matrix_row(&matrixZUY.matrix, i);
					gsl_blas_dcopy(&rowYTransS.vector, &rowZUY.vector);
					gsl_blas_dscal(gsl_matrix_get(&U.matrix, i, k),
							&rowZUY.vector);
				}

				/* Update Mu using BLAS */
				rowMu = gsl_matrix_row(&Mu.matrix, k);
				rowMu0 = gsl_matrix_row(&Mu0.matrix, k);
				gsl_vector_memcpy(&rowMu.vector, &rowMu0.vector);

				rowZ = gsl_matrix_column(&Z.matrix, k);
				gsl_blas_dgemv(CblasTrans, 1.0, &matrixYTransS.matrix,
						&rowZ.vector, *kappa0, &rowMu.vector);
				gsl_blas_dscal(1. / (gsl_vector_get(SumZU, k) + *kappa0),
						&rowMu.vector);
				/* Update Precision (cluster specific) */
				rowPrecision = gsl_matrix_row(&Precision.matrix, k);
				matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector,
						*py, *py);
				rowLambda0 = gsl_matrix_row(&Lambda0.matrix, k);
				matrixLambda0 = gsl_matrix_view_vector(&rowLambda0.vector, *py,
						*py);
				status2 = up_date_precision(&matrixZUY.matrix, &rowMu.vector,
						&matrixPrecision.matrix, gsl_vector_get(SumZ, k),
						gsl_vector_get(SumZU, k), DiagOne, &rowMu0.vector,
						&matrixLambda0.matrix, *kappa0, nu0[k]);
				/* Update Mixing Proportions */
				gsl_vector_set(&W.vector, k, gsl_vector_get(SumZ, k) / (*ly));
			}
		}

		/* Solve for nu */
		if ((*nuEstimate >= 1) && (iter < *BSolve)) // Only estimate nu in the first 100(BSolve) iterations
		{
			params2.SumZ = SumZ;
			params2.SumZU = SumZU;
			params2.SumZlogU = SumZlogU;
			params2.nu = nu;

			if (*nuEstimate == 1) {
				params2.K = -1;
				F2.params = &params2;

				/* Initialize the solver s2 to use the function F2 and the initial search interval */
				status = gsl_root_fsolver_set(s2, &F2, *nuLow, *nuUp);
				// Rprintf("K=%d  status=%d  %s\n",params2.K,status,gsl_strerror(status));
				if (status == 0) {
					iterSolve = 0;
					do {
						iterSolve++;
						/* Perform one iteration of Brent's algorithm */
						status = gsl_root_fsolver_iterate(s2);
						/* current estimate of nu */
						nuHat = gsl_root_fsolver_root(s2);
						/* current bracketing interval for solver */
						nu_lo = gsl_root_fsolver_x_lower(s2);
						nu_hi = gsl_root_fsolver_x_upper(s2);
						status = gsl_root_test_interval(nu_lo, nu_hi,
								*DiffSolve, 0);
					} while ((status == GSL_CONTINUE) && (iterSolve
							< *iterSolveMax));
					/* Set the new value of nu */
					for (k = 0; k < *K; k++)
						nu[k] = nuHat;
					// Rprintf("%d iterations, nu = %g\n",iterSolve,nuHat);
				}
			} else {
				for (k = 0; k < *K; k++) {
					params2.K = k;
					F2.params = &params2;

					/* Initialize the solver s2 to use the function F2 and the initial search interval */
					status = gsl_root_fsolver_set(s2, &F2, *nuLow, *nuUp);
					//Rprintf("K=%d  status=%d  %s\n",params2.K,status,gsl_strerror(status));
					if (status == 0) {
						iterSolve = 0;
						do {
							iterSolve++;
							/* Perform one iteration of Brent's algorithm */
							status = gsl_root_fsolver_iterate(s2);
							/* current estimate of nu */
							nuHat = gsl_root_fsolver_root(s2);
							/* current bracketing interval for solver */
							nu_lo = gsl_root_fsolver_x_lower(s2);
							nu_hi = gsl_root_fsolver_x_upper(s2);
							status = gsl_root_test_interval(nu_lo, nu_hi,
									*DiffSolve, 0);
						} while ((status == GSL_CONTINUE) && (iterSolve
								< *iterSolveMax));
						/* Set the new value of nu */
						nu[k] = nuHat;
						// Rprintf("%d iterations, nu = %g\n",iterSolve,nuHat);
					}
				}
			}
		} else if ((*nuEstimate == -1) && (iter < *BSolve)) // Only estimate nu in the first 100(BSolve) iterations
		{
			params3.W = &W.vector;
			params3.Mu = &Mu.matrix;
			params3.Precision = &Precision.matrix;
			params3.nu = nu;
			if (*transform > 0) {
				params3.lambda = lambda;
				if (*transform == 1)
					params3.YTrans = YTrans;
				else
					params3.YTrans = YTransS;
			}
			F3.params = &params3;

			/* Initialize the minimizer s3 to use the function F3 and the initial search interval */
			status = gsl_min_fminimizer_set(s3, &F3, *nu, *nuLow, *nuUp);
			if (status == 0) {
				iterSolve = 0;
				do {
					iterSolve++;
					/* Perform one iteration of Brent's algorithm */
					status = gsl_min_fminimizer_iterate(s3);
					/* current estimate of nu */
					nuHat = gsl_min_fminimizer_x_minimum(s3);
					/* current bracketing interval for minimizer */
					nu_lo = gsl_min_fminimizer_x_lower(s3);
					nu_hi = gsl_min_fminimizer_x_upper(s3);
					status = gsl_min_test_interval(nu_lo, nu_hi, *DiffSolve, 0);
				} while ((status == GSL_CONTINUE)
						&& (iterSolve < *iterSolveMax));
				/* Set the new value of nu */
				for (k = 0; k < *K; k++)
					nu[k] = nuHat;
			}
		}

		Diff = fabs(*logLike - logLikeOld) / fabs(logLikeOld);
		iter++;
	}

	/* One more E-step to compute the final z's and u's */
	if ((*transform <= 1) & (status2 == 0)) {
		up_date_z_u(logY, YTrans, &W.vector, &Mu.matrix, &Precision.matrix,
				&Z.matrix, &U.matrix, SumZ, SumZU, SumZlogU, nu, lambda,
				logLike, *transform, *nuEstimate, 1);
	} else if (status2 == 0) {
		up_date_z_uS(logY, YTransS, &W.vector, &Mu.matrix, &Precision.matrix,
				&Z.matrix, &U.matrix, SumZ, SumZU, SumZlogU, SumZlogY, nu,
				lambda, logLike, *nuEstimate, 1);
	}

	/** We encountered an error when computing the covariance matrix **/
	if (status2 != 0) {
		*logLike = GSL_NAN;
	}
	//Rprintf("End, logLikelihood=%f\n",iter,*logLike);

	/* Output Precision to be the covariance matrix */
	for (k = 0; k < *K; k++) {
		gsl_matrix_set_identity(DiagOne);
		rowPrecision = gsl_matrix_row(&Precision.matrix, k);
		matrixPrecision
		= gsl_matrix_view_vector(&rowPrecision.vector, *py, *py);
		gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0,
				&matrixPrecision.matrix, DiagOne);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0,
				&matrixPrecision.matrix);
	}

	/* Compute uncertainty and flagOutliers */
	if (*nuEstimate != 0) {
		for (k = 0; k < *K; k++)
			u_cutoff[k] = (nu[k] + *py) / (nu[k] + *py * gsl_cdf_fdist_Pinv(
					u_cutoff[k], *py, nu[k]));
	}
	for (i = 0; i < *ly; i++) {
		rowZ = gsl_matrix_row(&Z.matrix, i);
		label[i] = gsl_vector_max_index(&rowZ.vector);
		uncertainty[i] = 1 - gsl_vector_get(&rowZ.vector, label[i]);
		if (gsl_matrix_get(&U.matrix, i, label[i]) < u_cutoff[label[i]]
		                                                      || gsl_matrix_get(&Z.matrix, i, label[i]) < *z_cutoff)
			flagOutliers[i] = 1;
		//  else
		// flagOutliers[i]=0;
	}

	gsl_vector_free(SumZ);
	gsl_vector_free(SumZU);
	gsl_matrix_free(DiagOne);
	gsl_matrix_free(ZUY);
	gsl_matrix_free(YTrans);
	gsl_matrix_free(ZY);
	gsl_matrix_free(tmpOmega);
	gsl_matrix_free(MuGMatrix);

	if (*transform >= 1 || (*transform == 0 && *lambda != 1)) {
		gsl_matrix_free(logY);
		if (*transform >= 1)
			gsl_root_fsolver_free(s);
		if (*transform > 1) {
			gsl_matrix_free(YTransS);
			gsl_vector_free(SumZlogY);
		}
	}
	if(*transform==1&&*model==2){
		gsl_matrix_free(params.MuUpdate);
		gsl_matrix_free(params.PrecisionUpdate);
	}
	if (*nuEstimate >= 1) {
		gsl_vector_free(SumZlogU);
		gsl_root_fsolver_free(s2);
	} else if (*nuEstimate == -1)
		gsl_min_fminimizer_free(s3);
}

/* Compute the precision matrix and its cholesky decomposition */
int up_date_precision(gsl_matrix *ZUY, gsl_vector *Mu, gsl_matrix *Precision,
		double SumZ, double SumZU, gsl_matrix *DiagOne, gsl_vector *Mu0,
		gsl_matrix *Lambda0, double kappa0, double nu0) {
	int status = 0;
	int p = Mu->size;
	gsl_vector *muTmp;

	muTmp = gsl_vector_calloc(p);

	gsl_matrix_set_identity(DiagOne);
	gsl_matrix_memcpy(Precision, Lambda0);
	gsl_vector_memcpy(muTmp, Mu);
	gsl_vector_sub(muTmp, Mu0);

	/* Compute Prec=kappa0*(mu_g-mu0)%*%(mu_g-mu0)'+Prec */
	gsl_blas_dsyr(CblasLower, kappa0, muTmp, Precision);

	/* Compute Prec=-SumZU*mu_g%*%mu_g'+Prec */
	gsl_blas_dsyr(CblasLower, -SumZU, Mu, Precision);

	/* Compute Prec= sumzu yy' */
	gsl_blas_dsyrk(CblasLower, CblasTrans, 1. / (SumZ + nu0 + p + 1.), ZUY, 1.
			/ (SumZ + nu0 + p + 1.), Precision);
	/* Compute the cholesky decomposition of the covariance matrix = LL'*/
	status = gsl_linalg_cholesky_decomp(Precision);
	if (status != 0) {
		return (status);
	}
	/* Compute L'^{-1} */
	gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0,
			Precision, DiagOne);
	/* Compute Precision=L'^{-1}L^{-1} */
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0,
			Precision);
	/* Compute the cholesky decomposition of the precision matrix */
	status = gsl_linalg_cholesky_decomp(Precision);
	gsl_vector_free(muTmp);
	if (status != 0) {
		return (status);
	}
	return (0);
}

/* Compute Z, U, (SumZlogU) and also SumZ, SumZU, logLike */
void up_date_z_u(gsl_matrix *logY, gsl_matrix *YTrans, gsl_vector *W,
		gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U,
		gsl_vector *SumZ, gsl_vector *SumZU, gsl_vector *SumZlogU, double *nu,
		double *lambda, double *logLike, int transform, int nuEstimate,
		int last) {

	int i = 0, j = 0, k = 0, ly = (*YTrans).size1, py = (*YTrans).size2, K =
			(*Mu).size1;
	gsl_vector_view rowYTrans, rowMu, rowPrecision, rowZ;
	gsl_matrix_view matrixPrecision;
	double logJacobian = 0;
	double normConstant = 0; // normalizing constant for Z
	double tmpLike = 0; // W * density function
	double like = 0; // Sum of (W * density function) wrt one observation

	/*	Initialize elements to zero*/
	*logLike = 0;
	gsl_vector_set_zero(SumZ);
	gsl_vector_set_zero(SumZU);
	if (nuEstimate > 0)
		gsl_vector_set_zero(SumZlogU);

	for (i = 0; i < ly; i++) {
		if (transform == 1 || (transform == 0 && *lambda != 1)) {
			for (j = 0; j < py; j++)
				logJacobian += (*lambda - 1) * gsl_matrix_get(logY, i, j);
		}
		rowYTrans = gsl_matrix_row(YTrans, i);
		normConstant = 0;
		like = 0;
		for (k = 0; k < K; k++) {
			rowMu = gsl_matrix_row(Mu, k);
			rowPrecision = gsl_matrix_row(Precision, k);
			matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector, py,
					py);
			/* Here I assume the cholesky decomposition has been done */
			/* Compute W * density function */
			tmpLike = gsl_vector_get(W, k) * gsl_ran_mvnt_pdf(
					&rowYTrans.vector, &rowMu.vector, &matrixPrecision.matrix,
					nu[k], 1, 0);
			/* E[Z|y] (un-normalized) */
			gsl_matrix_set(Z, i, k, tmpLike);
			/* Compute the normalizing constant (for Z) */
			normConstant += gsl_matrix_get(Z, i, k);
			/* Compute E[u|y,z] */
			gsl_matrix_set(U, i, k, (py + nu[k]) / (gsl_pow_2(gsl_mahalanobis(
					&matrixPrecision.matrix, &rowYTrans.vector, &rowMu.vector,
					1)) + nu[k]));
			/*  Compute the likelihood wrt to one observation*/
			like += tmpLike;
		}
		*logLike += log(like);
		/* Scale the Z's so that they sum to one */
		rowZ = gsl_matrix_row(Z, i);
		gsl_blas_dscal(1.0 / normConstant, &rowZ.vector);

		for (k = 0; k < K; k++) {
			if (last == 0) // last=0 means not the final EM iteration
			{
				/* Compute E[ZlogU|y] (unadjusted) */
				if (nuEstimate > 0)
					gsl_vector_set(SumZlogU, k, gsl_vector_get(SumZlogU, k)
							+ gsl_matrix_get(Z, i, k) * gsl_sf_log(
									gsl_matrix_get(U, i, k)));

				/* Here I assume the cholesky decomposition has been done */
				/* Note that Z is used to store ZU */
				/* Note that U is used to store (ZU)^(1/2) */
				/* I store the square root for efficiency */
				gsl_vector_set(SumZ, k, gsl_vector_get(SumZ, k)
						+ gsl_matrix_get(Z, i, k));
				/* Compute ZU */
				gsl_matrix_set(Z, i, k, gsl_matrix_get(Z, i, k)
						* gsl_matrix_get(U, i, k));
				gsl_vector_set(SumZU, k, gsl_vector_get(SumZU, k)
						+ gsl_matrix_get(Z, i, k));
				/* Compute (ZU)^(1/2) */
				gsl_matrix_set(U, i, k, sqrt(gsl_matrix_get(Z, i, k)));
			}
		}
	}
	if ((nuEstimate > 0) && (last == 0)) { /* add adjustment term to E[ZlogU|y] */
		for (k = 0; k < K; k++)
			gsl_vector_set(SumZlogU, k, gsl_vector_get(SumZlogU, k)
					+ gsl_vector_get(SumZ, k) * (gsl_sf_psi((nu[k] + py) / 2)
							- gsl_sf_log((nu[k] + py) / 2)));
	}
	if (transform == 1 || (transform == 0 && *lambda != 1)) {
		*logLike += logJacobian;
	}
}

/* Compute Z, U, (SumZlogU, SumZlogY) and also SumZ, SumZU, logLike */
void up_date_z_uS(gsl_matrix *logY, gsl_matrix *YTransS, gsl_vector *W,
		gsl_matrix *Mu, gsl_matrix *Precision, gsl_matrix *Z, gsl_matrix *U,
		gsl_vector *SumZ, gsl_vector *SumZU, gsl_vector *SumZlogU,
		gsl_vector *SumZlogY, double *nu, double *lambda, double *logLike,
		int nuEstimate, int last) {
	int i = 0, j = 0, k = 0, ly = (*logY).size1, py = (*logY).size2, K =
			(*Mu).size1;
	gsl_vector_view rowYTransS, subRowYTransS, rowMu, rowPrecision, rowZ;
	gsl_matrix_view matrixPrecision;
	double logJacobian = 0;
	double normConstant = 0; // normalizing constant for Z
	double tmpLike = 0; // W * density function
	double like = 0; // Sum of (W * density function) wrt one observation

	/*	Initialize elements to zero*/
	*logLike = 0;
	gsl_vector_set_zero(SumZ);
	gsl_vector_set_zero(SumZU);
	if (nuEstimate > 0)
		gsl_vector_set_zero(SumZlogU);
	gsl_vector_set_zero(SumZlogY);

	for (i = 0; i < ly; i++) {
		rowYTransS = gsl_matrix_row(YTransS, i);
		normConstant = 0;
		like = 0;
		for (k = 0; k < K; k++) {
			logJacobian = 0;
			for (j = 0; j < py; j++)
				logJacobian += gsl_matrix_get(logY, i, j);
			logJacobian *= lambda[k] - 1;

			rowMu = gsl_matrix_row(Mu, k);
			rowPrecision = gsl_matrix_row(Precision, k);
			matrixPrecision = gsl_matrix_view_vector(&rowPrecision.vector, py,
					py);
			/* Here I assume the cholesky decomposition has been done */
			/* Compute W * density function */
			subRowYTransS
			= gsl_vector_subvector(&rowYTransS.vector, k * py, py);
			tmpLike = gsl_vector_get(W, k) * gsl_ran_mvnt_pdf(
					&subRowYTransS.vector, &rowMu.vector,
					&matrixPrecision.matrix, nu[k], 1, 0) * gsl_sf_exp(
							logJacobian);
			/* E[Z|y] (un-normalized) */
			gsl_matrix_set(Z, i, k, tmpLike);
			/* Compute the normalizing constant (for Z) */
			normConstant += gsl_matrix_get(Z, i, k);
			/* Compute E[u|y,z] */
			gsl_matrix_set(U, i, k, (py + nu[k]) / (gsl_pow_2(gsl_mahalanobis(
					&matrixPrecision.matrix, &subRowYTransS.vector,
					&rowMu.vector, 1)) + nu[k]));
			/*  Compute the likelihood wrt to one observation*/
			like += tmpLike;
		}
		*logLike += log(like);
		/* Scale the Z's so that they sum to one */
		rowZ = gsl_matrix_row(Z, i);
		gsl_blas_dscal(1.0 / normConstant, &rowZ.vector);

		for (k = 0; k < K; k++) {
			if (last == 0) // last=0 means not the final EM iteration
			{
				/* Compute Sum(ZlogY) */
				for (j = 0; j < py; j++)
					gsl_vector_set(SumZlogY, k, gsl_vector_get(SumZlogY, k)
							+ gsl_matrix_get(Z, i, k) * gsl_matrix_get(logY, i,
									j));
				/* Compute E[ZlogU|y] (unadjusted) */
				if (nuEstimate > 0)
					gsl_vector_set(SumZlogU, k, gsl_vector_get(SumZlogU, k)
							+ gsl_matrix_get(Z, i, k) * gsl_sf_log(
									gsl_matrix_get(U, i, k)));

				/* Here I assume the cholesky decomposition has been done */
				/* Note that Z is used to store ZU */
				/* Note that U is used to store (ZU)^(1/2) */
				/* I store the square root for efficiency */
				gsl_vector_set(SumZ, k, gsl_vector_get(SumZ, k)
						+ gsl_matrix_get(Z, i, k));
				/* Compute ZU */
				gsl_matrix_set(Z, i, k, gsl_matrix_get(Z, i, k)
						* gsl_matrix_get(U, i, k));
				gsl_vector_set(SumZU, k, gsl_vector_get(SumZU, k)
						+ gsl_matrix_get(Z, i, k));
				/* Compute (ZU)^(1/2) */
				gsl_matrix_set(U, i, k, sqrt(gsl_matrix_get(Z, i, k)));
			}
		}
	}
	if ((nuEstimate > 0) && (last == 0)) { /* add adjustment term to E[ZlogU|y] */
		for (k = 0; k < K; k++)
			gsl_vector_set(SumZlogU, k, gsl_vector_get(SumZlogU, k)
					+ gsl_vector_get(SumZ, k) * (gsl_sf_psi((nu[k] + py) / 2)
							- gsl_sf_log((nu[k] + py) / 2)));
	}
}

void getEstimates(double *y, int *ly, int *py, int *K, double *mu,
		double *precision, double *nu, double *z, double *u) {
	int i = 0, k = 0;
	gsl_vector *SumZ = gsl_vector_calloc(*K), *SumZU = gsl_vector_calloc(*K);
	gsl_matrix *Weight = gsl_matrix_calloc(*ly, *K); // Weight = Z*U
	gsl_matrix *TempMatrix = gsl_matrix_calloc(*py, *py);
	gsl_matrix *DiagOne = gsl_matrix_calloc(*py, *py);
	gsl_matrix *WeightedY = gsl_matrix_calloc(*ly, *py); // WeightedY = ZU*Y
	gsl_matrix_view Y, Mu, Precision, Z, U, matrixPrecision;
	gsl_vector_view rowMu, rowPrecision, rowWeightedY;

	/* Create matrix and vector views*/
	Y = gsl_matrix_view_array(y, *ly, *py);
	Mu = gsl_matrix_view_array(mu, *K, *py);
	Precision = gsl_matrix_view_array(precision, *K, *py * *py);
	Z = gsl_matrix_view_array(z, *ly, *K);
	U = gsl_matrix_view_array(u, *ly, *K);
	gsl_matrix_set_identity(DiagOne);

	/* Compute Weight = Z*U */
	gsl_matrix_memcpy(Weight, &U.matrix);
	gsl_matrix_mul_elements(Weight, &Z.matrix);
	for (i = 0; i < *ly; i++) {
		for (k = 0; k < *K; k++) {
			gsl_vector_set(SumZ, k, gsl_vector_get(SumZ, k) + gsl_matrix_get(
					&Z.matrix, i, k));
			gsl_vector_set(SumZU, k, gsl_vector_get(SumZU, k) + gsl_matrix_get(
					Weight, i, k));
		}
	}

	/* Compute Mu using BLAS */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Weight, &Y.matrix, 0.0,
			&Mu.matrix);
	for (k = 0; k < *K; k++) {
		rowMu = gsl_matrix_row(&Mu.matrix, k);
		gsl_blas_dscal(1. / gsl_vector_get(SumZU, k), &rowMu.vector);

		/* Compute WeightedY = ZU*Y */
		gsl_matrix_memcpy(WeightedY, &Y.matrix);
		for (i = 0; i < *ly; i++) {
			rowWeightedY = gsl_matrix_row(WeightedY, i);
			gsl_blas_dscal(gsl_matrix_get(Weight, i, k), &rowWeightedY.vector);
		}

		/* Compute Precision (cluster specific) */
		rowPrecision = gsl_matrix_row(&Precision.matrix, k);
		matrixPrecision
		= gsl_matrix_view_vector(&rowPrecision.vector, *py, *py);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0 / gsl_vector_get(SumZ, k),
				WeightedY, &Y.matrix, 0.0, TempMatrix);
		gsl_blas_dsyr(CblasLower, -gsl_vector_get(SumZU, k) / gsl_vector_get(
				SumZ, k), &rowMu.vector, TempMatrix);
		gsl_blas_dsymm(CblasLeft, CblasLower, nu[k] / (nu[k] - 2.0),
				TempMatrix, DiagOne, 0.0, &matrixPrecision.matrix);
	}

	gsl_vector_free(SumZ);
	gsl_vector_free(SumZU);
	gsl_matrix_free(Weight);
	gsl_matrix_free(TempMatrix);
	gsl_matrix_free(DiagOne);
	gsl_matrix_free(WeightedY);
}

// =====================================================
// = Update SigmaG, the MAP estimate of the covariance =
// =====================================================
int ECMUpdateSigmaG2(gsl_vector* U, gsl_matrix* YTrans,gsl_matrix* MuGMatrix, gsl_matrix* Precision,double SumZ,double nu0){
	gsl_matrix* sqrtZU = gsl_matrix_calloc(U->size,YTrans->size2);
	int py = YTrans->size2;
	gsl_matrix* DiagOne = gsl_matrix_calloc(py,py);
	gsl_matrix_set_identity(DiagOne);
	int i=0;
	int status=0;
	for (int i=0;i<YTrans->size2;i++){
		gsl_matrix_set_col(sqrtZU,i,U);
	};
	gsl_matrix_mul_elements(MuGMatrix,sqrtZU);
	gsl_matrix_mul_elements(sqrtZU,YTrans);
	gsl_matrix_sub(sqrtZU,MuGMatrix);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,sqrtZU,sqrtZU,1.0,Precision);
	gsl_matrix_scale(Precision,1.0/(SumZ + py + 1 + nu0));

	/*Cholesky decomposition of covariance */
	status = gsl_linalg_cholesky_decomp(Precision);
	if(status!=0){
		//printf("Error computing cholesky decomposition of covariance in  ECMUpdatesigmag2\n");
		//gsl_matrix_fprintf(stdout,Precision,"%f");
		return(status);
	}

	/* Compute L'^{-1} */
	status=gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0,
			Precision, DiagOne);

	/* Compute Precision=L'^{-1}L^{-1} */
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, Precision);

	/* Compute the cholesky decomposition of the precision matrix */
	status = gsl_linalg_cholesky_decomp(Precision);

	gsl_matrix_free(sqrtZU);
	gsl_matrix_free(DiagOne);
	return(0);
}

int ECMUpdateMUg(gsl_matrix* Precision, gsl_vector* Mu, gsl_vector* ZY,
		gsl_matrix *Omega0, double SumZU, gsl_matrix *DiagOne) {
	/* Precision variable has been initialized in flowClust.R  to be the inverse of the Cholesky decomposition of the precision matrix */
	int status = 0;
	int p = Mu->size;
	gsl_vector *tmpMu = gsl_vector_calloc(p);
	gsl_matrix_set_identity(DiagOne);
	gsl_matrix* tmpOmega = gsl_matrix_calloc(Omega0->size1,Omega0->size2);
	gsl_matrix* tmpPrecision = gsl_matrix_calloc(Precision->size1,Precision->size2);
	gsl_matrix_memcpy(tmpOmega,Omega0);
	/* Compute the precision matrix */
	gsl_matrix_memcpy(tmpPrecision, Precision);
	/* Invert the cholesky decomposition of the precision, square, decompose, invert, and square again */
	gsl_blas_dtrsm(CblasLeft, CblasUpper,CblasNoTrans, CblasNonUnit,1.0, tmpPrecision, DiagOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, tmpPrecision);
	gsl_linalg_cholesky_decomp(tmpPrecision);
	gsl_matrix_set_identity(DiagOne);
	gsl_blas_dtrsm(CblasLeft, CblasUpper,CblasNoTrans, CblasNonUnit,1.0, tmpPrecision, DiagOne);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DiagOne, DiagOne, 0.0, tmpPrecision);



	/* Compute Mu = Precision^{-1}*ZY + Mu */

	gsl_blas_dgemv(CblasNoTrans, 1.0, tmpPrecision, ZY, 1.0, Mu);

	/*Compute Omega0 = sumzu*matrixPrecision^{-1} + Omega0 */
	gsl_matrix_set_identity(DiagOne);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, SumZU, tmpPrecision, DiagOne, 1.0,
			tmpOmega);

	/*Invert tmpOmega) */
	gsl_linalg_cholesky_decomp(tmpOmega);
	status = gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
			1.0,tmpOmega,DiagOne);
	if(status!=0){
		Rprintf("Error inverting tmpOmega in MuUpdate\n");

	}
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,DiagOne,DiagOne,0.0,tmpOmega);
	gsl_matrix_set_identity(DiagOne);


	/* Finally rowMu = tmpOmega*rowMu */
	gsl_vector_memcpy(tmpMu, Mu);
	gsl_blas_dgemv(CblasNoTrans, 1.0, tmpOmega, tmpMu, 0.0, Mu);

	gsl_vector_free(tmpMu);
	gsl_matrix_free(tmpOmega);
	gsl_matrix_free(tmpPrecision);
	if (status != 0) {
		return (status);
	}

}

// =========================================================
// = Order the priors by their distance to each component. =
// =========================================================
void reorderPriors(gsl_matrix *Mu, gsl_matrix *Mu0, gsl_matrix *Lambda0, gsl_matrix *Omega0,double *nu0,double* w0, int* oorder) {
	int K = Mu->size1; //number of clusters in Mu
	int py = Mu->size2; //number of dimensions in Mu
	gsl_matrix* edist = gsl_matrix_calloc(K,K);

	int *order = calloc(K,sizeof(int));
	double *renu0 = calloc(K,sizeof(double));
	double *rew0 = calloc(K,sizeof(double));

	gsl_matrix *diff = gsl_matrix_calloc(K,py);
	gsl_vector_view rowMu0, rowMu, rowDiff;
	gsl_matrix *reMu0 = gsl_matrix_calloc(K,py);
	gsl_matrix *reLambda0=gsl_matrix_calloc(K,py*py);
	gsl_matrix *reOmega0 = gsl_matrix_calloc(K,py*py);
	int *reoorder = calloc(K,sizeof(int));

	for (int j=0;j < K;j++){
		for (int i=0;i<K;i++){
			rowMu0 = gsl_matrix_row(Mu0,j);
			rowMu = gsl_matrix_row(Mu,i);
			rowDiff = gsl_matrix_row(diff,i);
			gsl_vector_memcpy(&rowDiff.vector,&rowMu.vector);
			gsl_vector_sub(&rowDiff.vector,&rowMu0.vector);
			gsl_matrix_set(edist,j,i,gsl_blas_dnrm2(&rowDiff.vector));
		}
		//get the minimum element;
		gsl_vector_view edistView = gsl_matrix_row(edist,j);
		if(j==0){
			order[j] =  gsl_vector_min_index(&edistView.vector);
		}else{
			for(int i=0;i<j;i++){
				gsl_vector_set(&edistView.vector,order[i],GSL_POSINF);
			}
			order[j] = gsl_vector_min_index(&edistView.vector);
		}
	}	
	//Reorder Mu0 and Lambda0, nu0, and w0 according to order
	for(int i=0;i<K;i++){
		gsl_vector_view viewReMu0 = gsl_matrix_row(reMu0,order[i]);
		gsl_matrix_get_row(&viewReMu0.vector,Mu0,i); //copy row of Mu0 into the appropriate row of reMu0

		gsl_vector_view viewReLambda0 = gsl_matrix_row(reLambda0,order[i]);
		gsl_matrix_get_row(&viewReLambda0.vector,Lambda0,i);

		gsl_vector_view viewReOmega0 = gsl_matrix_row(reOmega0,order[i]);
		gsl_matrix_get_row(&viewReOmega0.vector,Omega0,i);
		renu0[order[i]]=nu0[i];
		reoorder[order[i]]=oorder[i];
		rew0[order[i]]=w0[i];
	}

	//copy reordered varaibles to the original variables.
	gsl_matrix_memcpy(Mu0,reMu0);
	gsl_matrix_memcpy(Lambda0,reLambda0);
	gsl_matrix_memcpy(Omega0,reOmega0);
	memcpy(nu0,renu0,K*sizeof(double));
	memcpy(oorder,reoorder,K*sizeof(int));
	memcpy(w0,rew0,K*sizeof(double));

	gsl_matrix_free(diff);
	gsl_matrix_free(reMu0);
	gsl_matrix_free(reLambda0);
	gsl_matrix_free(edist);
	gsl_matrix_free(reOmega0);
	free(order);
	free(rew0);
	free(renu0);
}

