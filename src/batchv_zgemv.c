/**
 * @file batchv_zgemv.c
 *
 * Part of API test for Batched BLAS routines.
 *
 * @author Samuel D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-06-01
 *
 * @precisions normal z -> c d s
 *
 **/

#include <cblas.h>
#include "bblas.h"

#define COMPLEX
void batchv_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int batch_count, int* info)
{
	/* Local variables */
//	int first_index = 0;
	int batch_iter = 0;
	char func_name[15] = "batchv_zgemv";

	if (batch_count < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
	}

	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	{
		/* Check input arguments */
		if ((trans[batch_iter] != BblasTrans) &&
			(trans[batch_iter] != BblasNoTrans) &&
			(trans[batch_iter] != BblasConjTrans))
		{
			xerbla_batch(func_name, BBLAS_ERR_TRANS, batch_iter);
			info[batch_iter] = BBLAS_ERR_TRANS;
		}

		if (m[batch_iter] < 0)
		{
			xerbla_batch(func_name, BBLAS_ERR_M, batch_iter);
			info[batch_iter] = BBLAS_ERR_M;
		}

		if (n[batch_iter] < 0)
		{
			xerbla_batch(func_name, BBLAS_ERR_N, batch_iter);
			info[batch_iter] = BBLAS_ERR_N;
		}

		/* Column major */
		if ((lda[batch_iter] < 1) && (lda[batch_iter] < m[batch_iter]))
		{
			xerbla_batch(func_name, BBLAS_ERR_LDA, batch_iter);
			info[batch_iter] = BBLAS_ERR_LDA;
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

		/* Call CBLAS */
		cblas_zgemv(
			BblasColMajor,
			trans[batch_iter],
			m[batch_iter], n[batch_iter],
			CBLAS_SADDR( alpha[batch_iter] ),
			arrayA[batch_iter], lda[batch_iter],
			arrayx[batch_iter], incx[batch_iter],
			CBLAS_SADDR( beta[batch_iter] ),
			arrayy[batch_iter], incy[batch_iter]);
		/* Successful */
		info[batch_iter] = BBLAS_SUCCESS;
	} /* End variable size for loop */
}
#undef COMPLEX
