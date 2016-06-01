/**
 * @file bblas_z.h
 *
 * @brief BBLAS header file for double _Complex routines.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#ifndef BBLAS_Z_H
#define BBLAS_Z_H


#include "bblas_types.h"
#include "bblas_macros.h"


#define COMPLEX

/*
 * Declarations of level 3 BATCH BLAS  - alphabetical order
 */


/*
 * Declarations of level 2 BATCH BLAS - alphabetical order
 */
void batch_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int batch_count, enum BBLAS_OPTS batch_opts,
	int* info);

void batchf_zgemv(
	const enum BBLAS_TRANS trans,
	const int m, const int n,
	const BBLAS_Complex64_t alpha,
	const BBLAS_Complex64_t **arrayA, const int lda,
	const BBLAS_Complex64_t **arrayx, const int incx,
	const BBLAS_Complex64_t beta,
	BBLAS_Complex64_t **arrayy, const int incy,
	const int batch_count, int* info)
	);

void batchg_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int group_count, const int *group_size,
	int* info);

void batchv_zgemv(
	const enum BBLAS_TRANS *trans,
	const int *m, const int *n,
	const BBLAS_Complex64_t *alpha,
	const BBLAS_Complex64_t **arrayA, const int *lda,
	const BBLAS_Complex64_t **arrayx, const int *incx,
	const BBLAS_Complex64_t *beta,
	BBLAS_Complex64_t **arrayy, const int *incy,
	const int batch_count, int* info);

/*
 * Error handler
 */
int xerbla_batch(char *func_name, int error, int subproblem_ind);

#undef COMPLEX
#endif /* BBLAS_Z_H */
