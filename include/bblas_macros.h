/**
 * @file bblas_macros.h
 *
 * @brief BBLAS macro definitions.
 *
 *  BBLAS is a software package provided by Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-02-20
 *
 * Contains macros for success and error codes, max, min, and passing complex values to CBLAS.
 *
 **/
#ifndef BBLAS_MACROS_H
#define BBLAS_MACROS_H


/*
 * BBLAS Return Codes
 */
#define BBLAS_SUCCESS                0 //!< Computation completed successfully
/* Batch error codes */
#define BBLAS_ERR_BATCH_COUNT       -1 //!< Error in batch_count
#define BBLAS_ERR_BATCH_OPTS        -2 //!< Error in batch_opts
/* Subproblem error codes */
#define BBLAS_ERR_M                 -3 //!< Error in M
#define BBLAS_ERR_N                 -4 //!< Error in N
#define BBLAS_ERR_K                 -5 //!< Error in K
#define BBLAS_ERR_LDA               -6 //!< Error in lda
#define BBLAS_ERR_LDB               -7 //!< Error in ldb
#define BBLAS_ERR_LDC               -8 //!< Error in ldc
#define BBLAS_ERR_UPLO              -9 //!< Error in uplo
#define BBLAS_ERR_TRANSA            -10 //!< Error in transA
#define BBLAS_ERR_TRANSB            -11 //!< Error in transB
#define BBLAS_ERR_TRANS             -12 //!< Error in trans
#define BBLAS_ERR_SIDE              -13 //!< Error in side
#define BBLAS_ERR_DIAG              -14 //!< Error in diag
#define BBLAS_ERR_INCX              -15 //!< Error in incx
#define BBLAS_ERR_INCY              -16 //!< Error in incy
#define BBLAS_ERR_STRIDEA           -17 //!< Error in strideA
#define BBLAS_ERR_STRIDEB           -17 //!< Error in strideB
#define BBLAS_ERR_STRIDEC           -17 //!< Error in strideC
/* Group error codes */
#define BBLAS_ERR_GROUP_COUNT       -30 //!<Error in group_count
#define BBLAS_ERR_GROUP_SIZE        -31 //!<Error in group_size

#define NAME_LENGTH 30 //!< Maximum length of routine name e.g. zgemm_batch

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b)) //!< Take max of two numbers
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b)) //!< Take min of two numbers
#endif

/**
 * CBLAS requires for complex scalar arguments to be passed by address rather than by value
 **/
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif

#endif /* BBLAS_MACROS_H */
