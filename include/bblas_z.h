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
void batch_zgemm(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB,
    const int *ldb, const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
    enum BBLAS_OPTS batch_opts, int *info);

void batchf_zgemm(
    const enum BBLAS_TRANS transA, const enum BBLAS_TRANS transB,
    const int M,  const int N, const int K,
    const BBLAS_Complex64_t alpha,
    const BBLAS_Complex64_t **arrayA, const int lda,
    const BBLAS_Complex64_t **arrayB,
    const int ldb, const BBLAS_Complex64_t beta,
    BBLAS_Complex64_t **arrayC, const int ldc,
    const int batch_count, int *info);

void batchg_zgemm(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc,
	const int group_count, const int *group_size,
    int *info);

void batchv_zgemm(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc,
    const int batch_count, int *info);

// One pointer (stride) approach

void batchf_zgemm_stride(
    const enum BBLAS_TRANS transA, const enum BBLAS_TRANS transB,
    const int M,  const int N, const int K,
    const BBLAS_Complex64_t alpha,
    const BBLAS_Complex64_t *arrayA, const int lda, const int strideA,
    const BBLAS_Complex64_t *arrayB, const int ldb, const int strideB,
    const BBLAS_Complex64_t beta,
    BBLAS_Complex64_t *arrayC, const int ldc, const int strideC,
    const int batch_count, int *info);

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
	const int batch_count, int info);

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
 * Declarations of level 1 BATCH BLAS - alphabetical order
 */
void batch_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const * x,
		  const int *incx,
		  BBLAS_Complex64_t const * const * y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, const enum BBLAS_OPTS batch_opts,
		  int * info);

void batchf_zdotu_sub(
		  const int n,
		  BBLAS_Complex64_t const * const * x,
		  const int incx,
		  BBLAS_Complex64_t const * const * y,
		  const int incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, int * info);

void batchg_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const * x,
		  const int *incx,
		  BBLAS_Complex64_t const * const * y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int group_count, const int *group_size,
		  int* info);

void batchv_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const * x,
		  const int *incx,
		  BBLAS_Complex64_t const * const * y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int batch_count, int* info);

#undef COMPLEX
#endif /* BBLAS_Z_H */
