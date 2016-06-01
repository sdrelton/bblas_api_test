/**
 *
 * @file bblas_types.h
 *
 * @brief BBLAS typedefs and enumerates.
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
 * Contains enumerates for parameter values and typedefs to support various platforms.
 *
 **/
#ifndef BBLAS_TYPES_H
#define BBLAS_TYPES_H

/*
 *  BBLAS enumerates for parameter values
 */

/**
 * Enum for row major or column major.
 **/
enum BBLAS_ORDER {BblasRowMajor=101, BblasColMajor=102};

/**
 * Enum for trans.
 **/
enum BBLAS_TRANS {BblasNoTrans=111, BblasTrans=112, BblasConjTrans=113};

/**
 * Enum for uplo.
 **/
enum BBLAS_UPLO {BblasUpper=121, BblasLower=122};

/**
 * Enum for diag.
 **/
enum BBLAS_DIAG {BblasNonUnit=131, BblasUnit=132};

/**
 * Enum for side.
 **/
enum BBLAS_SIDE {BblasLeft=141, BblasRight=142};

/**
 * Enum for batch_opts
 **/
enum BBLAS_OPTS{BBLAS_FIXED=0, BBLAS_VARIABLE= 1};

/*
 *  BBLAS typedefs
 */
typedef int  BBLAS_bool; //!< Define BBLAS_bool
typedef long BBLAS_index; //!< Define BBLAS_index
typedef long BBLAS_size; //!< Define BBLAS_size
typedef double BBLAS_Double_t; //!< Define BBLAS_Double_t


/*
 * BBLAS Complex numbers
 */


#define BBLAS_HAS_COMPLEX_H 1 //!< Is complex arithmetic is available?

#if defined(_WIN32)
# include <float.h>
# if defined(__INTEL_COMPILER)
    /* Fix name conflict within the cabs prototype (_Complex) that
     * conflicts with a C99 keyword.  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex
    #undef complex
typedef float  _Complex BBLAS_Complex32_t; //!< Define BBLAS_Complex32_t
    typedef double _Complex BBLAS_Complex64_t; //!< Define BBLAS_Complex64_t
# else
    /* Use MS VC complex class */
    #include <complex>
    typedef std::complex<float> BBLAS_Complex32_t; //!< Define BBLAS_Complex32_t
    typedef std::complex<double> BBLAS_Complex64_t; //!< Define BBLAS_Complex64_t
    #undef BBLAS_HAS_COMPLEX_H
# endif
# define isnan _isnan
# define isinf !_finite

#else /* defined(_WIN32) */

    typedef float  _Complex BBLAS_Complex32_t; //!< Define BBLAS_Complex32_t
    typedef double _Complex BBLAS_Complex64_t; //!< Define BBLAS_Complex64_t

#endif /* defined(_WIN32) */

/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef BBLAS_HAS_COMPLEX_H
#endif

#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define BBLAS_DEPRECATED  __attribute__((__deprecated__)) //!< Define BBLAS_DEPRECIATED
#else
#define BBLAS_DEPRECATED //!< Define BBLAS_DEPRECIATED
#endif /* __GNUC__ */

#ifdef BBLAS_HAS_COMPLEX_H
#include <complex.h>

#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because the names in C++ are name-mangled. */
#if !defined(_WIN32)
extern double cabs(BBLAS_Complex64_t z);
extern BBLAS_Complex64_t conj(BBLAS_Complex64_t z);
#endif
extern float cabsf(BBLAS_Complex32_t z);
extern double cimag(BBLAS_Complex64_t z);
extern double creal(BBLAS_Complex64_t z);

#ifdef __cplusplus
}
#endif
#endif /* defined(BBLAS_HAS_COMPLEX_H) */

#endif  /* BBLAS_TYPES_H */
