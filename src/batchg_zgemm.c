/**
 * @file batchg_zgemm.c
 *
 *  @brief BBLAS batchg_zgemm double _Complex routine.
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
    <b>batchg_zgemm</b> is a batch version of zgemm for groups 
    of batched with different sizes.
	It performs the matrix-matrix operations

	arrayC[i] = alpha[i]*op( arrayA[i] )*op( arrayB[i] ) + beta[i]*arrayC[i],

	where op( X ) is one of

	op( X ) = X      or
	op( X ) = X**T   or
	op( X ) = X**H,

	alpha[i] and beta[i] are scalars, and arrayA[i], arrayB[i] and C are matrices, with
	op( arrayA[i] ) an m by k matrix, op( arrayB[i] ) a k by n matrix and arrayC[i]
	an m by n matrix.

	Group Batch Operations
	----------------------
	
    - all parameters that are arrays must have length at least group_count.
	- all parameters that are arrays must have all values set.

	Parameters
	----------
	@param[in]
	transA   Array of <tt>enum BBLAS_TRANS</tt>.
	         On entry, transA[i] specifies the form of op( arrayA[i] ) to be used in
			 the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayA[i] ) = arrayA[i].
      -     = BblasTrans:      op( arrayA[i] ) = arrayA[i]**T.
      -     = BblasConjTrans:  op( arrayA[i] ) = arrayA[i]**H.

	@param[in]
	transB   Array of <tt>enum BBLAS_TRANS</tt>.
	         On entry, transB[i] specifies the form of op( arrayB[i] ) to be used in
             the matrix multiplication as follows:
      -     = BblasNoTrans:    op( arrayB[i] ) = arrayB[i].
      -     = BblasTrans:      op( arrayB[i] ) = arrayB[i]**T.
      -     = BblasConjTrans:  op( arrayB[i] ) = arrayB[i]**H.

	@param[in]
	M       Array of <tt>int</tt>.
            Each element M[i] specifies the number of rows of the matrix
            op( arrayA[i] ) and of the matrix arrayC[i]. M[i] must be greater than zero.

	@param[in]
	N       Array of <tt>int</tt>.
            Each element N[i] specifies the number of columns of the matrix
            op( arrayB[i] ) and the number of columns of the matrix arrayC[i].
 		    N[i] must be greater than zero.

    @param[in]
    K      Array of <tt>int</tt>.
	       Each element K[i] specifies the number of columns of the matrix
           op( arrayA[i] ) and the number of rows of the matrix op( arrayB[i] ).
		   K[i] must be greater than zero.

    @param[in]
    alpha   Array of <tt>complex_16</tt>.

    @param[in]
    arrayA  Array of pointers.
            Each element arrayA[i] is a pointer to a COMPLEX_16 matrix of
 		    dimension lda[i] by Ka[i], where Ka[i] is K[i] when transA[i] = BblasNoTrans,
 		    and is M[i] otherwise.
            When using transA[i] = BblasNoTrans the leading M[i] by K[i]
            part of arrayA[i] must contain the matrix elements, otherwise
            the leading  K[i] by M[i] part of arrayA[i] must contain the
            matrix elements.

    @param[in]
    lda     Array of <tt>int</tt>.
            Each element lda[i] specifies the first dimension of arrayA[i] as declared
            in the calling (sub) program. When transA[i] = BblasNoTrans then
            lda[i] must be at least max( 1, M[i] ), otherwise lda[i] must be at
            least max( 1, K[i] ).

    @param[in]
    arrayB  Array of pointers.
  		    Each element arrayB[i] is a pointer to a COMPLEX_16 matrix of
  		    dimension ldb[i] by Kb[i], where Kb[i] is N[i] when transB[i] = BblasNoTrans,
  		    and is K[i] otherwise.
            When using transB[i] = BblasNoTrans the leading K[i] by N[i]
            part of arrayB[i] must contain the matrix elements, otherwise
            the leading N[i] by K[i] part of arrayB[i] must contain the
            matrix elements.


    @param[in]
    ldb     Array of <tt>int</tt>.
  		    Each element ldb[i] specifies the first dimension of arrayB[i] as declared
            in the calling (sub) program. When transB[i] = BblasNoTrans then
            ldb[i] must be at least max( 1, K[i] ), otherwise ldb[i] must be at
            least max( 1, N[i] ).

    @param[in]
    beta    Array of <tt>complex_16</tt>.
            When beta[i] is set to zero arrayC[i] need not be set on input.

    @param[in,out]
    arrayC  Array of pointers.
	        Each element arrayC[i] is a pointer to a COMPLEX_16 matrix of
      	    dimension ldc[i] by N[i].
  	        Before entry, the leading M[i] by N[i] part of the arrayC[i] must
            contain a matrix C, except when beta is zero, in which
            case C need not be set on entry.
            On exit, the matrix arrayC[i] is overwritten by the M[i] by N[i] matrix
            ( alpha[i]*op( arrayA[i] )*op( arrayB[i] ) + beta[i]*arrayC[i] ).

    @param[in]
    ldc     Array of <tt>int</tt>.
  	        Each element ldc[i] specifies the first dimension of arrayC[i] as declared
  	        in the calling (sub) program. The value ldc[i] must be at least
  	        max( 1, M[i] )

    @param[in]
    batch_count  <tt>int</tt>
                 The number of matrices to operate on.

    @param[in,out]
    info    Array of <tt>int</tt>.
            Each element info[i] is the error return code of the ith group,
            these need not be set on entry.
  	    The error codes can be found in bblas_macros.h.
 **/
void batchg_zgemm(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex64_t *alpha,
    const BBLAS_Complex64_t **arrayA, const int *lda,
    const BBLAS_Complex64_t **arrayB,
    const int *ldb, const BBLAS_Complex64_t *beta,
    BBLAS_Complex64_t **arrayC, const int *ldc, 
	const int group_count, const int *group_size,
    int *info)
{
    /* Local variables */
    int group_iter;
    int offset = 0;
    char func_name[15] = "batchg_zgemm";

	/* Check input arguments */
    if (group_count < 0)
    {
	    xerbla_batch(func_name, BBLAS_ERR_BATCH_COUNT, -1);
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

        /* Call to bblas_zgemm_batchf */
        batchf_zgemm(
            transA[group_iter],
            transB[group_iter],
            M[group_iter],
            N[group_iter],
            K[group_iter],
            alpha[group_iter],
            arrayA+offset,
            lda[group_iter],
            arrayB+offset,
            ldb[group_iter],
            beta[group_iter],
            arrayC+offset,
            ldc[group_iter],
            group_size[group_iter],
            &info[group_iter]);    

        offset += group_size[group_iter];    
    }
}
#undef COMPLEX
