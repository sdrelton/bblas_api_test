/**
 * @file batchf_zdotu_sub.c
 *
 * Part of API test for Batched BLAS routines.
 *
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date 
 *
 * @precisions normal z -> c 
 *
 **/

#include <cblas.h>
#include "bblas.h"

#define COMPLEX

void batchf_zdotu_sub(
		  const int n,
		  BBLAS_Complex64_t const * const * x,
		  const int incx,
		  BBLAS_Complex64_t const * const * y,
		  const int incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, int *info)
{
	/* Local variables */
	int first_index = 0;
	int batch_iter = 0;
	char func_name[15] = "batchf_zdotu";

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

	if (n < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_N, first_index);
		info[first_index] = BBLAS_ERR_N;

	}


	if (incx < 1)
	{
		xerbla_batch(func_name, BBLAS_ERR_INCX, first_index);
		info[first_index] = BBLAS_ERR_INCX;
	}
	if (incy < 1)
	{
		xerbla_batch(func_name, BBLAS_ERR_INCY, first_index);
		info[first_index] = BBLAS_ERR_INCY;
	}
	/* Call CBLAS */
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	{
	  cblas_zdotu_sub(
		      n,
		      (void *)x[batch_iter],
		      incx,
		      (void *)y[batch_iter],
		      incy,
		      &dotu[batch_iter]);
	  /* Successful */
	} /* End fixed size for loop */
	info[first_index] = BBLAS_SUCCESS;
}
#undef COMPLEX
