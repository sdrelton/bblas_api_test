#include "example_z.h"
#include <stdlib.h>
#include <lapacke.h>

void random_zvec(int len, BBLAS_Complex64_t *vec)
{
//	BBLAS_Complex64_t *work =
//		(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*len);
//	int info = 0;
	double scalar = 1.0;
	int seed = 2;
	int zero = 0;
	int one = 1;
	int lenm1 = len-1;
	int colmaj = BblasColMajor;
	LAPACKE_zlagge(colmaj, len, one, lenm1, zero, &scalar, vec, len, &seed);
}

void random_mat(int m, int n, BBLAS_Complex64_t *mat)
{
//	BBLAS_Complex64_t *work = malloc(sizeof(BBLAS_Complex64_t)*m*n);
//	int info = 0;
	double scalar = 1.0;
	int seed = 2;
	int mm1 = m-1;
	int nn1 = n-1;
	int colmaj = BblasColMajor;
	LAPACKE_zlagge(colmaj, m, n, mm1, nn1, &scalar, mat, m, &seed);
}

BBLAS_Complex64_t randz()
{
	BBLAS_Complex64_t val;
	random_mat(1, 1, &val);
	return val;
}

void set_params_fixed_zgemv(struct zgemv_batchf_example *zgemv_example)
{
	int batch_count = rand() % 50;
	zgemv_example->batch_count = batch_count;
	zgemv_example->trans = rand() % 3+111;
	int m = rand() % 200;
	zgemv_example->m = m;
	int n = rand() % 200;
	zgemv_example->n = n;
	zgemv_example->alpha = randz();
	/* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*n);
		random_mat(m, n, arrayA[i]);
	}
	zgemv_example->arrayA = arrayA;
	zgemv_example->lda = zgemv_example->m;
	/* arrayx */
	BBLAS_Complex64_t **arrayx =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	if (zgemv_example->trans == BblasNoTrans)
	{
		curdim = n;
	}
	else
	{
		curdim = m;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayx[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayx[i]);
	}
	zgemv_example->arrayx = arrayx;
	zgemv_example->incx = 1;
	zgemv_example->beta = randz();
	/* arrayy */
	BBLAS_Complex64_t **arrayy =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	if (zgemv_example->trans == BblasNoTrans)
	{
		curdim = m;
	}
	else
	{
		curdim = n;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayy[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayy[i]);
	}
	zgemv_example->arrayy = arrayy;
	zgemv_example->incy = 1;
}

void set_params_variable_zgemv(struct zgemv_batchv_example *zgemv_example)
{
	int batch_count = rand() % 50;
	zgemv_example->batch_count = batch_count;
	int *trans = (int*) malloc(sizeof(int)*batch_count);
	for(int i = 0; i < batch_count; i++)
	{
		trans[i] = rand() % 3 + 111;
	}
	zgemv_example->trans = (enum BBLAS_TRANS*) trans;
	int *m = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		m[i] = rand() % 200;
	}
	zgemv_example->m = m;
	int *n = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		n[i] = rand() % 200;
	}
	zgemv_example->n = n;
	BBLAS_Complex64_t *alpha = malloc(sizeof(BBLAS_Complex64_t)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		alpha[i] = randz();
	}
	zgemv_example->alpha = alpha;
	/* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m[i]*n[i]);
		random_mat(m[i], n[i], arrayA[i]);
	}
	zgemv_example->arrayA = arrayA;
	int *lda = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		lda[i] = m[i];
	}
	zgemv_example->lda = lda;
	/* arrayx */
	BBLAS_Complex64_t **arrayx =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	for (int i = 0; i < batch_count; i++)
	{
		if (trans[i] == BblasNoTrans)
		{
			curdim = n[i];
		}
		else
		{
			curdim = m[i];
		}
		arrayx[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayx[i]);
	}
	zgemv_example->arrayx = arrayx;
	int *incx = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		incx[i] = 1;
	}
	zgemv_example->incx = incx;
	BBLAS_Complex64_t *beta = malloc(sizeof(BBLAS_Complex64_t)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		beta[i] = randz();
	}
	zgemv_example->beta = beta;
	/* arrayy */
	BBLAS_Complex64_t **arrayy =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		if (trans[i] == BblasNoTrans)
		{
			curdim = m[i];
		}
		else
		{
			curdim = n[i];
		}
		arrayy[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayy[i]);
	}
	zgemv_example->arrayy = arrayy;
	int *incy = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		incy[i] = 1;
	}
	zgemv_example->incy = incy;
}
