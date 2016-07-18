/**
 * @file batchv_zdotu_sub.c
 *
 * Part of API test for Batched BLAS routines.
 *
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date 
 *
 * @precisions normal z -> c d s
 *
 **/

#include <cblas.h>
#include "bblas.h"

#define COMPLEX

void batchv_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const * x,
		  const int *incx,
		  BBLAS_Complex64_t const * const * y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, int* info)
{
	/* Local variables */
	int batch_iter = 0;
	char func_name[15] = "batchv_zdotu";

	/*initialize the result */
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	  {
	    dotu[batch_iter] = (BBLAS_Complex64_t)0.0;
	  }
	
	/* Check input arguments */
	if (batch_count < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
	}

	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	  {
	    if (n[batch_iter] < 0)
	      {
		xerbla_batch(func_name, BBLAS_ERR_N, batch_iter);
		info[batch_iter] = BBLAS_ERR_N;
	      }
	    if (incx[batch_iter] < 1)
	      {
		xerbla_batch(func_name, BBLAS_ERR_INCX, batch_iter);
		info[batch_iter] = BBLAS_ERR_INCX;
	      }
	    if (incy[batch_iter] < 1)
	      {
		xerbla_batch(func_name, BBLAS_ERR_INCY, batch_iter);
		info[batch_iter] = BBLAS_ERR_INCY;
	      }

	    cblas_zdotu_sub(
		      n[batch_iter],
		      (void *)x[batch_iter],
		      incx[batch_iter],
		      (void *)y[batch_iter],
		      incy[batch_iter],
		      &dotu[batch_iter]);
	    /* Successful */
	    info[batch_iter] = BBLAS_SUCCESS;
	  } /* End fixed size for loop */
}
#undef COMPLEX
