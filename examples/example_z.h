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

void random_zvec(int len, BBLAS_Complex64_t *vec);

void random_mat(int m, int n, BBLAS_Complex64_t *mat);

void set_params_fixed_zgemv(struct zgemv_batchf_example *zgemv_example);
