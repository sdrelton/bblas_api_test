#include "example_z.h"
#include <math.h>
#include <stdlib.h>
#include <lapacke.h>

void random_zvec(int len, BBLAS_Complex64_t *vec)
{
	BBLAS_Complex64_t* work = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*len);
	int info = 0;
	lapacke_zlagge(len, 1, len-1, 0, 1.0, vec, len, 2, work, info);
	free(work);
}

void random_mat(int m, int n, BBLAS_Complex64_t *mat)
{
	BBLAS_Complex64_t* work = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*n);
	int info = 0;
	lapacke_zlagge(n, n, m-1, n-1, 1.0, mat, m, 2, work, info);
	free(work);
}

BBLAS_Complex64_t randz()
{
	BBLAS_Complex64_t val;
	random_mat(1, 1, &val);
	return val;
}

void set_params_fixed_zgemv(struct zgemv_batchf_example *zgemv_example)
{
	int batch_count = rand() % 200;
	zgemv_example->batch_count = batch_count;
	zgemv_example->trans = rand() % 3;
	int m = rand() % 200;
	zgemv_example->m = m;
	int n = rand() % 200;
	zgemv_example->n = n;
	zgemv_example->alpha = randz();
	/* arrayA */
	BBLAS_Complex64_t **arrayA = (BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*n);
		random_mat(m, n, arrayA[i]);
	}
	zgemv_example->arrayA = arrayA;
	zgemv_example->lda = zgemv_example->m;
	/* arrayx */
	BBLAS_Complex64_t **arrayx = (BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	if (trans == BblasNoTrans)
	{
		curdim = n;
	}
	else
	{
		curdim = m;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayx[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_vec(curdim, arrayx[i]);
	}
	zgemv_example->arrayx = arrayx;
	zgemv_example->incx = 1;
	zgemv_example->beta = randz();
		/* arrayx */
	BBLAS_Complex64_t **arrayy = (BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	if (trans == BblasNoTrans)
	{
		curdim = m;
	}
	else
	{
		curdim = n;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayy[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_vec(curdim, arrayy[i]);
	}
	zgemv_batchf_example->arrayy = arrayy;
	zgemv_batchf_example->incy = 1;


}
