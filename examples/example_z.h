#include "bblas_z.h"

struct zgemv_batchf_example {
	enum BBLAS_TRANS trans;
	int m;  int n;
	BBLAS_Complex64_t alpha;
	BBLAS_Complex64_t **arrayA;  int lda;
	BBLAS_Complex64_t **arrayx;  int incx;
	BBLAS_Complex64_t beta;
	BBLAS_Complex64_t **arrayy;  int incy;
	int batch_count;
};

struct zgemv_batchv_example {
	enum BBLAS_TRANS *trans;
	int *m; int *n;
	BBLAS_Complex64_t *alpha;
	BBLAS_Complex64_t **arrayA; int *lda;
	BBLAS_Complex64_t **arrayx; int *incx;
	BBLAS_Complex64_t *beta;
	BBLAS_Complex64_t **arrayy; int *incy;
	int batch_count;
};

struct zgemm_batchf_example {
    enum BBLAS_TRANS transA; 
    enum BBLAS_TRANS transB;
    int m; int n; int k;
    BBLAS_Complex64_t alpha;
    BBLAS_Complex64_t **arrayA; int lda;
    BBLAS_Complex64_t **arrayB; int ldb;  
    BBLAS_Complex64_t beta;
    BBLAS_Complex64_t **arrayC; int ldc; 
    int batch_count;
};

// Struct for the one pointer (stride) approach

struct zgemm_batchf_example_stride {
	enum BBLAS_TRANS transA; 
    enum BBLAS_TRANS transB;
    int m; int n; int k;
    BBLAS_Complex64_t alpha;
    BBLAS_Complex64_t *arrayA; int lda; int strideA;
    BBLAS_Complex64_t *arrayB; int ldb; int strideB; 
    BBLAS_Complex64_t beta;
    BBLAS_Complex64_t *arrayC; int ldc; int strideC; 
    int batch_count;
};

struct zgemm_batchv_example {
    enum BBLAS_TRANS *transA; 
    enum BBLAS_TRANS *transB;
    int *m; int *n; int *k;
    BBLAS_Complex64_t *alpha;
    BBLAS_Complex64_t **arrayA; int *lda;
    BBLAS_Complex64_t **arrayB; int *ldb;  
    BBLAS_Complex64_t *beta;
    BBLAS_Complex64_t **arrayC;  int *ldc; 
    int batch_count; 
};

void random_zvec(int len, BBLAS_Complex64_t *vec);

void random_mat(int m, int n, BBLAS_Complex64_t *mat);

void set_params_fixed_zgemv(struct zgemv_batchf_example *zgemv_example);

void set_params_variable_zgemv(struct zgemv_batchv_example *zgemv_example);

void set_params_fixed_zgemm(struct zgemm_batchf_example *zgemm_example);

void set_params_fixed_zgemm_stride(struct zgemm_batchf_example_stride *zgemm_example);

void set_params_variable_zgemm(struct zgemm_batchv_example *zgemm_example);
