/*
 * BoxCoxPriors.c
 *
 *  Created on: Oct 21, 2010
 *      Author: finak
 */
#include <string.h>
#include "flowClust.h"

double BoxCoxGradientPriors(double x, void *params) {
	struct BoxCox_params *p = (struct BoxCox_params *) params;
	gsl_matrix *Y = p->Y;
	gsl_vector *W = p->W; // will be updated
	gsl_matrix *Mu = p->Mu; // will be updated
	gsl_matrix *Precision = p->Precision; // will be updated
	gsl_matrix *Z = p->Z;
	gsl_matrix *U = p->U;
	gsl_matrix *logY = p->logY;
	gsl_matrix *YTrans = p->YTrans; // will be updated
	gsl_matrix *ZUY = p->ZUY; // will be updated
	gsl_vector *SumZ = p->SumZ;
	gsl_vector *SumZU = p->SumZU;
	gsl_vector *SumZlogY = p->SumZlogY;
	gsl_matrix *DiagOne = p->DiagOne;
	gsl_matrix *Mu0 = p->Mu0;
	gsl_matrix *Omega0 = p->Omega0;
	gsl_matrix *Lambda0 = p->Lambda0;
	int model = *p->model;
	int *oorder = p->oorder;
	double *w0 = p->w0;
	int iter = p->iter;
	int *itersolve = p->itersolve;
	double *nu0 = p->nu0, kappa0 = p->kappa0;

	int K1 = p->K;

	int i = 0, ly = Y->size1;
	int j = 0, py = Mu->size2;
	int k = 0, K = Mu->size1;
	int K0 = 0;
	
	double Yoriginal = 0, Tmp = 0, Tmp1 = 0, Tmp2 = 0, Tmp3 = 0;
	double logJacobian = 0; // derivative of log(Jacobian) wrt lambda
	double logLike = 0; // derivative of loglikelihood wrt lambda
	gsl_matrix *Y1 = gsl_matrix_alloc(ly, py); // derivative of YTrans wrt lambda
	gsl_matrix *ZUY1 = gsl_matrix_alloc(ly, K * py);
	gsl_matrix *PrecisionYTrans = gsl_matrix_alloc(ly, py);
	gsl_vector *Y1Sum = gsl_vector_alloc(py);
	gsl_matrix *tmpOmega = gsl_matrix_calloc(Omega0->size1,Omega0->size2);
	gsl_matrix *tmpLambda0 = gsl_matrix_calloc(K,py*py);
	gsl_matrix *tmpMu0 = gsl_matrix_calloc(K,py);
	gsl_matrix *MuGMatrix = gsl_matrix_calloc(ly, py);
	gsl_matrix *ZY = gsl_matrix_calloc(K, py); // ZU*YTrans
	gsl_vector *SMu = gsl_vector_alloc(py);
	gsl_matrix_view matrixZUY, matrixZUY1, matrixPrecision, matrixLambda0;
	gsl_vector_view rowZUY, rowZUY1, rowPrecisionYTrans, rowMu, rowMu0,
	rowPrecision, rowZ, rowLambda0;
	double w0sum=0;
	
	for(int ii=0;ii<K1;ii++){
		w0sum=w0sum+w0[ii];
	}
	if(K1!=-1){
		//Rprintf("K1=%d",K1);
		error("Estimation of lambda with priors is not implemented for cluster specific degrees of freedom or cluster specific lambda.\n");
	}
	for (i = 0; i < ly; i++) {
		for (j = 0; j < py; j++) {
			/* Compute YTrans */
			Yoriginal = gsl_matrix_get(Y, i, j);
			Tmp1 = sgn(Yoriginal) * pow(fabs(Yoriginal), x);
			gsl_matrix_set(YTrans, i, j, (Tmp1 - 1.0) / x);
			/* Compute Y1 */
			Tmp2 = gsl_matrix_get(logY, i, j) * x - 1.0;
			Tmp3 = (Tmp1 * Tmp2 + 1.0) / gsl_pow_2(x);
			gsl_matrix_set(Y1, i, j, Tmp3);
			logJacobian += gsl_matrix_get(logY, i, j);
		}
	}

	/* Update ZUY-matrix and compute ZUY1 */
	for (k = K0; k < K; k++) {
		matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * py, ly, py);
		gsl_matrix_memcpy(&matrixZUY.matrix, YTrans);
		matrixZUY1 = gsl_matrix_submatrix(ZUY1, 0, k * py, ly, py);
		gsl_matrix_memcpy(&matrixZUY1.matrix, Y1);

		for (i = 0; i < ly; i++) {
			/* ZUY = sqrt(ZU) times YTrans */
			rowZUY = gsl_matrix_row(&matrixZUY.matrix, i);
			gsl_blas_dscal(gsl_matrix_get(U, i, k), &rowZUY.vector);
			/* ZUY1 = (ZU) times Y1 */
			rowZUY1 = gsl_matrix_row(&matrixZUY1.matrix, i);
			gsl_blas_dscal(gsl_matrix_get(Z, i, k), &rowZUY1.vector);
		}
	}

	/* Update Mu  */
	gsl_matrix_memcpy(p->PrecisionUpdate,p->Precision);

	/* Use copies of the priors */
	gsl_matrix_memcpy(tmpOmega, Omega0);
	gsl_matrix_memcpy(tmpMu0,Mu0);
	gsl_matrix_memcpy(tmpLambda0,Lambda0);
	double* tmpnu0=calloc(K,sizeof(double));
	memcpy(tmpnu0,nu0,K*sizeof(double));
	int* tmpoorder=calloc(K,sizeof(int));
	double* tmpw0 = calloc(K,sizeof(double));
	memcpy(tmpoorder,oorder,K*sizeof(int));
	memcpy(tmpw0,w0,K*sizeof(double));
	
	/*Assume Omega0 is already inverted*/
	/*  Mu = Omega0^{-1} * Mu0*/
	//gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmpMu0, tmpOmega,
	//		0.0, p->MuUpdate);
	/* ZY = Z*Y */// dim(k*p)
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Z, YTrans, 0.0, ZY);
	for (k = K0; k < K; k++) {
		gsl_vector_view rowOmega0 = gsl_matrix_row(tmpOmega,k);
		gsl_matrix_view matrixOmega0 = gsl_matrix_view_vector(&rowOmega0.vector,py,py);
		rowMu = gsl_matrix_row(p->MuUpdate, k);
		rowMu0=gsl_matrix_row(tmpMu0,k);

		gsl_blas_dgemv(CblasNoTrans, 1.0,
				&matrixOmega0.matrix, &rowMu0.vector, 0.0, &rowMu.vector);


		rowPrecision = gsl_matrix_row(p->PrecisionUpdate, k);
		matrixPrecision
		= gsl_matrix_view_vector(&rowPrecision.vector, py, py);

		gsl_vector_view rowZY = gsl_matrix_row(ZY, k);
		ECMUpdateMUg(&matrixPrecision.matrix, &rowMu.vector, &rowZY.vector,
				&matrixOmega0.matrix, gsl_vector_get(SumZU, k), DiagOne);
	}
	reorderPriors(p->MuUpdate,tmpMu0,tmpLambda0,tmpOmega,tmpnu0,tmpw0,tmpoorder);
	/* Update Sigma */
	for (k = K0; k < K; k++) {
		rowPrecision = gsl_matrix_row(p->PrecisionUpdate, k);
		matrixPrecision
		= gsl_matrix_view_vector(&rowPrecision.vector, py, py);
		rowLambda0 = gsl_matrix_row(tmpLambda0, k);
		matrixLambda0 = gsl_matrix_view_vector(&rowLambda0.vector, py, py);
		matrixZUY = gsl_matrix_submatrix(ZUY, 0, k * py, ly, py);
		rowMu = gsl_matrix_row(p->MuUpdate, k);
		if (gsl_matrix_isnull(MuGMatrix)) {
			for (i = 0; i < ly; i++) {
				gsl_matrix_set_row(MuGMatrix, i, &rowMu.vector);
			}
		}
		gsl_matrix_memcpy(&matrixPrecision.matrix, &matrixLambda0.matrix);
		gsl_vector_view zu = gsl_matrix_column(U, k);
		ECMUpdateSigmaG2(&zu.vector, YTrans, MuGMatrix,
				&matrixPrecision.matrix, gsl_vector_get(SumZ, k), nu0[k]);
		gsl_matrix_set_zero(MuGMatrix);
	}
	for (k = K0; k < K; k++) {
		/* Update Mixing Proportions */
		gsl_vector_set(W, k, gsl_vector_get(SumZ, k) / ly);
		gsl_vector_set(W, k, (gsl_vector_get(SumZ, k)+w0[k])
				/ (w0sum+(ly)));

		rowMu = gsl_matrix_row(p->MuUpdate, k);

		/* Compute PrecisionYTrans = YTrans times sqrt(Precision) */
		gsl_matrix_memcpy(PrecisionYTrans, YTrans);
		gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0,
				&matrixPrecision.matrix, PrecisionYTrans);

		/* Compute Y1Sum = sqrt(Precision) times sum(ZUY1) */
		gsl_vector_set_zero(Y1Sum);
		for (i = 0; i < ly; i++) {
			rowZUY1 = gsl_vector_view_array(ZUY1->data + i * py * K + k * py,
					py);
			gsl_blas_daxpy(1.0, &rowZUY1.vector, Y1Sum);
		}
		gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit,
				&matrixPrecision.matrix, Y1Sum);

		/* Compute ZUY1 = ZUY1 times sqrt(Precision) */
		matrixZUY1 = gsl_matrix_submatrix(ZUY1, 0, k * py, ly, py);
		gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0,
				&matrixPrecision.matrix, &matrixZUY1.matrix);

		/* Compute PrecisionYTrans times ZUY1*/
		for (i = 0; i < ly; i++) {
			rowPrecisionYTrans = gsl_matrix_row(PrecisionYTrans, i);
			rowZUY1 = gsl_vector_view_array(ZUY1->data + i * py * K + k * py,
					py);
			gsl_blas_ddot(&rowPrecisionYTrans.vector, &rowZUY1.vector, &Tmp);
			logLike -= Tmp;
		}

		/* Compute SMu = Mu times Precision */
		gsl_blas_dcopy(&rowMu.vector, SMu);
		gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit,
				&matrixPrecision.matrix, SMu);
		/* Compute SMu times Y1Sum */
		gsl_blas_ddot(SMu, Y1Sum, &Tmp);
		logLike += Tmp;
	}

	logLike += logJacobian;
	//Rprintf("logLikeDeriv=%f\n",logLike);

	gsl_vector_free(Y1Sum);
	gsl_vector_free(SMu);
	gsl_matrix_free(ZY);
	gsl_matrix_free(Y1);
	gsl_matrix_free(tmpOmega);
	gsl_matrix_free(tmpLambda0);
	gsl_matrix_free(tmpMu0);
	gsl_matrix_free(PrecisionYTrans);
	gsl_matrix_free(ZUY1);
	gsl_matrix_free(MuGMatrix);
	free(tmpnu0);
	free(tmpoorder);
	free(tmpw0);
	return (logLike);
}
