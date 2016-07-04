/**
 * @file batchf_zgemm.c
 *
 *  @brief BBLAS batchf_zgemm double _Complex routine.
 *
 *  BBLAS is a software package provided by Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-06-03
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include <cblas.h>
#include "bblas.h"

#define COMPLEX

/**
    Purpose
    -------
    <b>batchf_zgemm</b> is a batch version of zgemm for fixed size.
	It performs the matrix-matrix operations

	arrayC[i] = alpha[i]*op( arrayA[i] )*op( arrayB[i] ) + beta[i]*arrayC[i],

	where op( X ) is one of

	op( X ) = X      or
	op( X ) = X**T   or
	op( X ) = X**H,

	alpha[i] and beta[i] are scalars, and arrayA[i], arrayB[i] and C are matrices, with
	op( arrayA[i] ) an m by k matrix, op( arrayB[i] ) a k by n matrix and arrayC[i]
	an m by n matrix.


	Fixed Batch Operations
	-----------------------------------

	- all parameters that are arrays must have length at least batch_count.
	- all parameters that are arrays must have all values set.

	the values of transA, transB, M, N, K, alpha, beta, lda, ldb, and ldc 
    are used for all computations.


	Parameters
	----------
	@param[in]
	transA   <tt>enum BBLAS_TRANS</tt>.
	         transA specifies the form of op( arrayA[i] ) to be used in
			 the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayA[i] ) = arrayA[i].
      -     = BblasTrans:      op( arrayA[i] ) = arrayA[i]**T.
      -     = BblasConjTrans:  op( arrayA[i] ) = arrayA[i]**H.

	@param[in]
	transB   <tt>enum BBLAS_TRANS</tt>.
	         On entry, transB specifies the form of op( arrayB[i] ) to be used in
             the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayB[i] ) = arrayB[i].
      -     = BblasTrans:      op( arrayB[i] ) = arrayB[i]**T.
      -     = BblasConjTrans:  op( arrayB[i] ) = arrayB[i]**H.

	@param[in]
	M       <tt>int</tt>.
            M specifies the number of rows of the matrix
            op( arrayA[i] ) and of the matrix arrayC[i]. 
            M must be greater than zero.

	@param[in]
	N       <tt>int</tt>.
            N specifies the number of columns of the matrix
            op( arrayB[i] ) and the number of columns of the matrix arrayC[i].
 		    N must be greater than zero.

    @param[in]
    K      <tt>int</tt>.
	       K specifies the number of columns of the matrix
           op( arrayA[i] ) and the number of rows of the matrix op( arrayB[i] ).
		   K must be greater than zero.

    @param[in]
    alpha   <tt>complex_16</tt>.

    @param[in]
    arrayA  Array of pointers.
            Each element arrayA[i] is a pointer to a COMPLEX_16 matrix of
 		    dimension lda by Ka, where Ka is K when transA = BblasNoTrans,
 		    and is M otherwise.
            When using transA = BblasNoTrans the leading M by K
            part of arrayA[i] must contain the matrix elements, otherwise
            the leading K by M part of arrayA[i] must contain the
            matrix elements.

    @param[in]
    lda     <tt>int</tt>.
            lda specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When transA = BblasNoTrans then
            lda must be at least max( 1, M ), otherwise lda must be at
            least max( 1, K ).

    @param[in]
    arrayB  Array of pointers.
  		    Each element arrayB[i] is a pointer to a COMPLEX_16 matrix of
  		    dimension ldb by Kb, where Kb is N when transB = BblasNoTrans,
  		    and is K otherwise.
            When using transB = BblasNoTrans the leading K by N
            part of arrayB[i] must contain the matrix elements, otherwise
            the leading N by K part of arrayB[i] must contain the
            matrix elements.


    @param[in]
    ldb     <tt>int</tt>.
  		    ldb specifies the first dimension of arrayB[i] as declared
            in the calling (sub) program. When transB = BblasNoTrans then
            ldb must be at least max( 1, K ), otherwise ldb must be at
            least max( 1, N ).

    @param[in]
    beta    <tt>complex_16</tt>.
            When beta is set to zero arrayC[i] need not be set on input.

    @param[in,out]
    arrayC  Array of pointers.
	        Each element arrayC[i] is a pointer to a COMPLEX_16 matrix of
      	    dimension ldc by N.
  	        Before entry, the leading M by N part of the arrayC[i] must
            contain a matrix C, except when beta is zero, in which
            case C need not be set on entry.
            On exit, the matrix arrayC[i] is overwritten by the M by N matrix
            ( alpha*op( arrayA[i] )*op( arrayB[i] ) + beta*arrayC[i] ).

    @param[in]
    ldc     <tt>int</tt>.
  	        Each element ldc specifies the first dimension of arrayC[i] as declared
  	        in the calling (sub) program. The value ldc must be at least
  	        max( 1, M )

    @param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

    @param[in,out]
    info    <tt>int</tt>.
            info is the error return code of the batch.
  	        The error codes can be found in bblas_macros.h.
 **/
void batchf_zgemm(
    const enum BBLAS_TRANS transA, const enum BBLAS_TRANS transB,
    const int M,  const int N, const int K,
    const BBLAS_Complex64_t alpha,
    const BBLAS_Complex64_t **arrayA, const int lda,
    const BBLAS_Complex64_t **arrayB,
    const int ldb, const BBLAS_Complex64_t beta,
    BBLAS_Complex64_t **arrayC, const int ldc, 
    const int batch_count, int *info)
{
    /* Local variables */
    int batch_iter;
    int LDA,  LDB;
    char func_name[15] = "batchf_zgemm";

	/* Check input arguments */
    if (batch_count < 0)
    {
	    xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
    }
    if ((transA != BblasNoTrans) &&
        (transA != BblasTrans) &&
        (transA != BblasConjTrans))
    {
        xerbla_batch(func_name, BBLAS_ERR_TRANSA, -1);
        *info  = BBLAS_ERR_TRANSA;
        return;
    }
    if ((transB != BblasNoTrans) &&
        (transB != BblasTrans) &&
        (transB != BblasConjTrans))
    {
        xerbla_batch(func_name, BBLAS_ERR_TRANSB, -1);
        *info = BBLAS_ERR_TRANSB;
        return;
    }
    if ( transA == BblasNoTrans )
    {
        LDA = M;
    } else
    {
        LDA = K;
    }
    if ( transB == BblasNoTrans )
    {
        LDB = K;
    } else
    {
        LDB = N;
    }
    if (M < 0)
    {
        xerbla_batch(func_name, BBLAS_ERR_M, -1);
        *info = BBLAS_ERR_M;
        return;
    }
    if (N < 0) {
        xerbla_batch(func_name, BBLAS_ERR_N, -1);
        *info = BBLAS_ERR_N;
        return;
    }
    if (K < 0) {
        xerbla_batch(func_name, BBLAS_ERR_K, -1);
        *info = BBLAS_ERR_K;
        return;
    }
    if (lda < max(1, LDA)) {
        xerbla_batch(func_name, BBLAS_ERR_LDA, -1);
        *info =  BBLAS_ERR_LDA;
        return;
    }
    if (ldb < max(1, LDB)) {
        xerbla_batch(func_name, BBLAS_ERR_LDB, -1);
        *info = BBLAS_ERR_LDB;
        return;
    }
    if (ldc < max(1, M)) {
        xerbla_batch(func_name, BBLAS_ERR_LDC, -1);
        *info = BBLAS_ERR_LDC;
        return;
    }
    /* Skip subproblems where nothing needs to be done */
    if (M == 0 || N == 0     ||
	    ((alpha == (BBLAS_Complex64_t)0.0 || K == 0) &&
	    beta == (BBLAS_Complex64_t)1.0 ))
    {
        *info =  BBLAS_SUCCESS;
        return;
    }
    for (batch_iter = 0; batch_iter < batch_count; batch_iter++)
    {
        /* Call to cblas_zgemm */
        cblas_zgemm(
        BblasColMajor,
        transA,
        transB,
        M,
        N,
        K,
        CBLAS_SADDR(alpha),
        arrayA[batch_iter],
        lda,
        arrayB[batch_iter],
        ldb,
        CBLAS_SADDR(beta),
        arrayC[batch_iter],
        ldc);
    }
    /* Successful */
    *info = BBLAS_SUCCESS;
}
#undef COMPLEX
