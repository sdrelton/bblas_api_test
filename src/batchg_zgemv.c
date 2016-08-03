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

void batchg_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int group_count, const int *group_size,
	int* info)
{
	/* Local variables */
	char func_name[15] = "batchg_zgemv";
	int group_iter = 0;
	int offset = 0; // How many subproblems solved so far

	/* Check group_count */
	if (group_count < 0)
	{
		xerbla_batch(func_name, BBLAS_ERR_GROUP_COUNT, -1);
		return;
	}

	/* Check group_size and call fixed batch computation */
	for (group_iter = 0; group_iter < group_count; group_iter++)
	{
		if (group_size[group_iter] < 0)
		{
			xerbla_batch(func_name, BBLAS_ERR_GROUP_SIZE, group_iter);
			info[group_iter] = BBLAS_ERR_GROUP_SIZE;
			continue;
		}

		/* Call fixed batch computation on the group */
		batchf_zgemv(
			trans[group_iter],
			m[group_iter], n[group_iter],
			alpha[group_iter],
			arrayA+offset, lda[group_iter],
			arrayx+offset, incx[group_iter],
			beta[group_iter],
			arrayy+offset, incy[group_iter],
			group_size[group_iter], info[group_iter]);

		/* Update offset */
		offset += group_size[group_iter];
	} // End of group loop
}
#undef COMPLEX
