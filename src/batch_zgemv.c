/**
 * @file batchg_zgemv.c
 *
 * Part of API test for Batched BLAS routines.
 *
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date 2016-06-01
 *
 * @precisions normal z -> c d s
 *
 **/

#include <cblas.h>
#include "bblas.h"

#define COMPLEX
void batch_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int batch_count, const enum BBLAS_OPTS batch_opts,
	int* info)
{
	/* Local variables */
	char func_name[15] = "batch_zgemv";

	/* Check input arguments */
	if (batch_count < 0)
    {
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
    }

	if (batch_opts == BBLAS_FIXED)
	{
		/* Call fixed size code */
		batchf_zgemv(trans[0],
					 m[0], n[0],
					 alpha[0],
					 arrayA, lda[0],
					 arrayx, incx[0],
					 beta[0],
					 arrayy, incy[0],
					 batch_count, info[0]);
	}
	else if (batch_opts == BBLAS_VARIABLE)
	  {
		/* Call variable size code */
		batchv_zgemv(trans,
					 m, n,
					 alpha,
					 arrayA, lda,
					 arrayx, incx,
					 beta,
					 arrayy, incy,
					 batch_count, info);
	}
	else
	{
		xerbla_batch(func_name, BBLAS_ERR_BATCH_OPTS, -1);
	}
}
