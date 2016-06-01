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
#include "auxiliary.h"


#define COMPLEX

/*
 *  Declarations of level 3 BATCH BLAS  - alphabetical order
 */
void bblas_zgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB,
    const int *ldb, const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_zhemm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_zsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void bblas_zsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_zher2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB, const int *ldb,
    const double  *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_zsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_zherk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const double  *beta, BBLAS_Complex64_t **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_ztrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    BBLAS_Complex64_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void bblas_ztrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    BBLAS_Complex64_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

int xerbla_batch(char *func_name, int error, int subproblem_ind);

#undef COMPLEX
#endif /* BBLAS_Z_H */
