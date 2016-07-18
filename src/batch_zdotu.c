/**
 * @file batch_zdotu_sub.c
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

void batch_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const *x,
		  const int *incx,
		  BBLAS_Complex64_t const * const *y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, const enum BBLAS_OPTS batch_opts,
		  int* info)
{
	/* Local variables */
	int first_index = 0;
	char func_name[15] = "batch_zdotu";

	/*initialize the result */
	for (int batch_iter = 0; batch_iter < batch_count; batch_iter++)
	  {
	    dotu[batch_iter] = (BBLAS_Complex64_t)0.0;
	  }
	
	/* Check input arguments */
	if (batch_count < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
	}

	if (batch_opts == BBLAS_FIXED)
	  {
	    /* Call fixed size code */
	    batchf_zdotu_sub(
			     n[first_index],
			     x,
			     incx[first_index],
			     y,
			     incy[first_index],
			     dotu,
			     batch_count, info);
	  }
	else if (batch_opts == BBLAS_VARIABLE)
	  {
	    /* Call variable size code */
	    batchv_zdotu_sub(
			     n,
			     x,
			     incx,
			     y,
			     incy,
			     dotu,
			     batch_count, info);
	    
	  }
	else
	  {
	    xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
	  }
}
#undef COMPLEX
