#include "flowClust.h"

double BoxCoxGradient(double x, void *params)
{
  struct BoxCox_params *p = (struct BoxCox_params *) params;     
  gsl_matrix *Y = p->Y;
  gsl_vector *W = p->W;      // will be updated
  gsl_matrix *Mu = p->Mu;    // will be updated
  gsl_matrix *Precision = p->Precision;  // will be updated
  gsl_matrix *Z = p->Z;
  gsl_matrix *U = p->U;
  gsl_matrix *logY = p->logY;
  gsl_matrix *YTrans = p->YTrans;    // will be updated
  gsl_matrix *ZUY = p->ZUY;   // will be updated
  gsl_vector *SumZ = p->SumZ;    
  gsl_vector *SumZU = p->SumZU;                
  gsl_vector *SumZlogY = p->SumZlogY;
  gsl_matrix *DiagOne = p->DiagOne;
  int K1 = p->K;

  int i=0, ly=Y->size1;    
  int j=0, py=Mu->size2;
  int k=0, K=Mu->size1;
  int K0=0;

  double Yoriginal=0,Tmp=0,Tmp1=0,Tmp2=0,Tmp3=0;
  double logJacobian=0;    // derivative of log(Jacobian) wrt lambda
  double logLike=0;   // derivative of loglikelihood wrt lambda
  gsl_matrix *Y1 = gsl_matrix_alloc(ly,py);    // derivative of YTrans wrt lambda
  gsl_matrix *ZUY1 = gsl_matrix_alloc(ly,K*py);
  gsl_matrix *PrecisionYTrans = gsl_matrix_alloc(ly,py);
  gsl_vector *Y1Sum = gsl_vector_alloc(py);     
  gsl_vector *SMu = gsl_vector_alloc(py);
  gsl_matrix_view matrixZUY, matrixZUY1, matrixPrecision; 
  gsl_vector_view rowZUY, rowZUY1, rowPrecisionYTrans, rowMu, rowPrecision, rowZ; 
  // Rprintf("K=%d  lambda=%lf\n",K1,x);

  for(i=0;i<ly;i++)
  {
    for(j=0;j<py;j++)
    {
            /* Compute YTrans */
      Yoriginal=gsl_matrix_get(Y,i,j);
      Tmp1=sgn(Yoriginal)*pow(fabs(Yoriginal),x);
      gsl_matrix_set(YTrans,i,j,(Tmp1-1.0)/x);
            /* Compute Y1 */
      Tmp2=gsl_matrix_get(logY,i,j)*x-1.0;
      Tmp3=(Tmp1*Tmp2+1.0)/gsl_pow_2(x);
      gsl_matrix_set(Y1,i,j,Tmp3);

      if (K1==-1)
        logJacobian+=gsl_matrix_get(logY,i,j);
    }
  }
  if (K1>-1)
    logJacobian=gsl_vector_get(SumZlogY,K1);

    /* Update ZUY-matrix and compute ZUY1 */
  if(K1>-1) K0=K1;
  for(k=K0;k<K;k++)
  {
    matrixZUY=gsl_matrix_submatrix(ZUY, 0, k*py, ly, py);
    gsl_matrix_memcpy(&matrixZUY.matrix, YTrans);
    matrixZUY1=gsl_matrix_submatrix(ZUY1, 0, k*py, ly, py);
    gsl_matrix_memcpy(&matrixZUY1.matrix, Y1);

    for(i=0;i<ly;i++)
    {
      /* ZUY = sqrt(ZU) times YTrans */
      rowZUY=gsl_matrix_row(&matrixZUY.matrix,i);
      gsl_blas_dscal(gsl_matrix_get(U,i,k), &rowZUY.vector);
      /* ZUY1 = (ZU) times Y1 */
      rowZUY1=gsl_matrix_row(&matrixZUY1.matrix,i);
      gsl_blas_dscal(gsl_matrix_get(Z,i,k), &rowZUY1.vector);
    }

    if(K1>-1) break;
  }

  /* Update Mu using BLAS */
  if(K1==-1)
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Z, YTrans, 0, Mu);
  for(k=K0;k<K;k++)
  { 
    rowMu=gsl_matrix_row(Mu,k);
    if(K1>-1)
    {
      rowZ=gsl_matrix_column(Z,k);
      gsl_blas_dgemv(CblasTrans, 1.0, YTrans, &rowZ.vector, 0, &rowMu.vector);
    }
    gsl_blas_dscal(1./gsl_vector_get(SumZU,k),&rowMu.vector);

    /* Update Precision (cluster specific) */
    rowPrecision=gsl_matrix_row(Precision,k);
    matrixPrecision=gsl_matrix_view_vector(&rowPrecision.vector,py, py);            
    matrixZUY=gsl_matrix_submatrix(ZUY, 0, k*py, ly, py);  
    up_date_precision(&matrixZUY.matrix, &rowMu.vector, &matrixPrecision.matrix, gsl_vector_get(SumZ,k), gsl_vector_get(SumZU,k), DiagOne);
    /* Update Mixing Proportions */
    gsl_vector_set(W,k,gsl_vector_get(SumZ,k)/ly);

    /* Compute PrecisionYTrans = YTrans times sqrt(Precision) */
    gsl_matrix_memcpy(PrecisionYTrans, YTrans);
    gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &matrixPrecision.matrix, PrecisionYTrans);                

    /* Compute Y1Sum = sqrt(Precision) times sum(ZUY1) */
    gsl_vector_set_zero(Y1Sum);
    for(i=0;i<ly;i++)
    {
      rowZUY1=gsl_vector_view_array(ZUY1->data+i*py*K+k*py,py);
      gsl_blas_daxpy(1.0, &rowZUY1.vector, Y1Sum);            
    }
    gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &matrixPrecision.matrix, Y1Sum);

    /* Compute ZUY1 = ZUY1 times sqrt(Precision) */
    matrixZUY1=gsl_matrix_submatrix(ZUY1, 0, k*py, ly, py);               
    gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, &matrixPrecision.matrix, &matrixZUY1.matrix);        

    /* Compute PrecisionYTrans times ZUY1*/
    for(i=0;i<ly;i++)
    {
      rowPrecisionYTrans=gsl_matrix_row(PrecisionYTrans,i);
      rowZUY1=gsl_vector_view_array(ZUY1->data+i*py*K+k*py,py);
      gsl_blas_ddot(&rowPrecisionYTrans.vector, &rowZUY1.vector, &Tmp);
      logLike-=Tmp;
    }

    /* Compute SMu = Mu times Precision */
    gsl_blas_dcopy(&rowMu.vector, SMu);
    gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &matrixPrecision.matrix, SMu);
    /* Compute SMu times Y1Sum */
    gsl_blas_ddot(SMu, Y1Sum, &Tmp);
    logLike+=Tmp;

    if(K1>-1) break;
  }

  logLike+=logJacobian;     

  gsl_vector_free(Y1Sum);
  gsl_vector_free(SMu);
  gsl_matrix_free(Y1);
  gsl_matrix_free(PrecisionYTrans);
  gsl_matrix_free(ZUY1);

    // Rprintf("logLike = %lf\n",logLike);
  return(logLike);
}


int sgn(double x)
{
  if (x >= 0) return 1; else return -1;
}
