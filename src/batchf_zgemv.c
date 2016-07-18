/**
 * @file batchf_zgemv.c
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

void batchf_zgemv(
	const enum BBLAS_TRANS trans,
	const int m, const int n,
	const BBLAS_Complex64_t alpha,
	const BBLAS_Complex64_t **arrayA, const int lda,
	const BBLAS_Complex64_t **arrayx, const int incx,
	const BBLAS_Complex64_t beta,
	BBLAS_Complex64_t **arrayy, const int incy,
	const int batch_count, int info)
{
	/* Local variables */
	int first_index = 0;
	int batch_iter = 0;
	char func_name[15] = "batchf_zgemv";

	/* Check input arguments */
	if (batch_count < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
	}

	if ((trans != BblasTrans) &&
		(trans != BblasNoTrans) &&
		(trans != BblasConjTrans))
	{
		xerbla_batch(func_name, BBLAS_ERR_TRANS, first_index);
		info = BBLAS_ERR_TRANS;
	}

	if (m < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_M, first_index);
		info = BBLAS_ERR_M;
	}

	if (n < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_N, first_index);
		info = BBLAS_ERR_N;
	}

	/* Column major */
	if ((lda < 1) && (lda < m))
	{
		xerbla_batch(func_name, BBLAS_ERR_LDA, first_index);
		info = BBLAS_ERR_LDA;
	}

	if (incx < 1)
	{
		xerbla_batch(func_name, BBLAS_ERR_INCX, first_index);
		info = BBLAS_ERR_INCX;
	}
	if (incy < 1)
	{
		xerbla_batch(func_name, BBLAS_ERR_INCY, first_index);
		info = BBLAS_ERR_INCY;
	}

	/* Call CBLAS */
	for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
	{
		cblas_zgemv(
			BblasColMajor,
			trans,
			m, n,
			CBLAS_SADDR( alpha ),
			arrayA[batch_iter], lda,
			arrayx[batch_iter], incx,
			CBLAS_SADDR( beta ),
			arrayy[batch_iter], incy);
	} /* End fixed size for loop */
	/* Successful */
	info = BBLAS_SUCCESS;
}
#undef COMPLEX
