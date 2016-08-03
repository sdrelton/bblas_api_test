/**
 *
 * @file test_gemv.c
 *
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-06-15
 * @precisions normal z -> c
 *
 **/
#include "bblas_z.h"
#include "example_z.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    printf("Testing gemv using some very simple tests.\n");

	/* Create simple test matrices */
	printf("Creating inputs for testing grouped zgemv...\n");

	/* Create two fixed size batches to use for fixed/grouped API */
	struct zgemv_batchf_example *zgemvf_example =
		(struct zgemv_batchf_example*) malloc(sizeof(struct zgemv_batchf_example));
	set_params_fixed_zgemv(zgemvf_example);

	struct zgemv_batchf_example *zgemvf_example2 =
		(struct zgemv_batchf_example*) malloc(sizeof(struct zgemv_batchf_example));
	set_params_fixed_zgemv(zgemvf_example2);

	/* Creat one variable size batch for variable API */
	struct zgemv_batchv_example *zgemvv_example =
		(struct zgemv_batchv_example*) malloc(sizeof(struct zgemv_batchv_example));
	set_params_variable_zgemv(zgemvv_example);

	printf("Finished creating inputs, running each API...\n");

	int infof = -1;
	int *infov = malloc(sizeof(int)*zgemvv_example->batch_count);

	printf("Computing with first API...\n");
	/*
	   First API: Use separate calls for fixed/variable

	   batchf_zgemv(
	   const BBLAS_TRANS trans,
	   const int m, const int n,
	   const Complex_ double alpha,
	   const Complex_ double **arrayA, const int lda,
	   const Complex_ double **arrayx, const int incx,
	   const Complex_ double beta,
	         Complex_ double **arrayy, const int incy,
	   const int batch_count, int *info)
	*/
	batchf_zgemv(
		(const enum BBLAS_TRANS) zgemvf_example->trans,
		(const int) zgemvf_example->m, (const int) zgemvf_example->n,
		(const BBLAS_Complex64_t) zgemvf_example->alpha,
		(const BBLAS_Complex64_t**) zgemvf_example->arrayA,
		(const int) zgemvf_example->lda,
		(const BBLAS_Complex64_t **) zgemvf_example->arrayx,
		(const int) zgemvf_example->incx,
		(const BBLAS_Complex64_t) zgemvf_example->beta,
		(BBLAS_Complex64_t **) zgemvf_example->arrayy,
		(const int) zgemvf_example->incy,
		(const int) zgemvf_example->batch_count, infof);

	/*
	  batchv_zgemv(
	  const BBLAS_TRANS *trans,
	  const int *m, const int *n,
	  const Complex_ double *alpha,
	  const Complex_ double **arrayA, const int *lda,
	  const Complex_ double **arrayx, const int *incx,
	  const Complex_ double *beta,
	        Complex_ double **arrayy, const int *incy,
	  const int batch_count, int *info)
	*/
	batchv_zgemv(
		(const enum BBLAS_TRANS*) zgemvv_example->trans,
		(const int*) zgemvv_example->m, (const int*) zgemvv_example->n,
		(const BBLAS_Complex64_t*) zgemvv_example->alpha,
		(const BBLAS_Complex64_t**) zgemvv_example->arrayA,
		(const int*) zgemvv_example->lda,
		(const BBLAS_Complex64_t **) zgemvv_example->arrayx,
		(const int*) zgemvv_example->incx,
		(const BBLAS_Complex64_t*) zgemvv_example->beta,
		(BBLAS_Complex64_t **) zgemvv_example->arrayy,
		(const int*) zgemvv_example->incy,
		(const int) zgemvv_example->batch_count, infov);

	printf("Computing with second API...\n");
	/* Group two fixed batched together */
	int group_trans[2] = {zgemvf_example->trans, zgemvf_example2->trans};
	int group_m[2] = {zgemvf_example->m, zgemvf_example2->m};
	int group_n[2] = {zgemvf_example->n, zgemvf_example2->n};
	BBLAS_Complex64_t group_alpha[2] = {zgemvf_example->alpha, zgemvf_example2->alpha};
	BBLAS_Complex64_t* group_arrayA[zgemvf_example->batch_count + zgemvf_example2->batch_count];
	for (int i = 0; i < zgemvf_example->batch_count; i++)
	{
		group_arrayA[i] = zgemvf_example->arrayA[i];
	}
	for (int i = 0; i < zgemvf_example2->batch_count; i++)
	{
		group_arrayA[i+zgemvf_example->batch_count] = zgemvf_example2->arrayA[i];
	}
	int group_lda[2] = {zgemvf_example->lda, zgemvf_example2->lda};
	BBLAS_Complex64_t* group_arrayx[zgemvf_example->batch_count + zgemvf_example2->batch_count];
	for (int i = 0; i < zgemvf_example->batch_count; i++)
	{
		group_arrayx[i] = zgemvf_example->arrayx[i];
	}
	for (int i = 0; i < zgemvf_example2->batch_count; i++)
	{
		group_arrayx[i+zgemvf_example->batch_count] = zgemvf_example2->arrayx[i];
	}
	int group_incx[2] = {zgemvf_example->incx, zgemvf_example2->incx};
	int group_beta[2] = {zgemvf_example->beta, zgemvf_example2->beta};
	BBLAS_Complex64_t* group_arrayy[zgemvf_example->batch_count + zgemvf_example2->batch_count];
	for (int i = 0; i < zgemvf_example->batch_count; i++)
	{
		group_arrayy[i] = zgemvf_example->arrayy[i];
	}
	for (int i = 0; i < zgemvf_example2->batch_count; i++)
	{
		group_arrayy[i+zgemvf_example->batch_count] = zgemvf_example2->arrayy[i];
	}
	int group_incy[2] = {zgemvf_example->incy, zgemvf_example2->incy};
	int group_count = 2;
	int group_size[2] = {zgemvf_example->batch_count, zgemvf_example2->batch_count};
	int infog[2];
	/*
	  batchg_zgemv(
	  const enum BBLAS_TRANS *trans,
	  const int *m, const int *n,
	  const BBLAS_Complex64_t *alpha,
	  const BBLAS_Complex64_t **arrayA, const int *lda,
	  const BBLAS_Complex64_t **arrayx, const int *incx,
	  const BBLAS_Complex64_t *beta,
	        BBLAS_Complex64_t **arrayy, const int *incy,
	  const int group_count, const int *group_size,
	  int* info)
	 */
	batchg_zgemv(
		(const enum BBLAS_TRANS*) group_trans,
		(const int*) group_m, (const int*) group_n,
		(const BBLAS_Complex64_t*) group_alpha,
		(const BBLAS_Complex64_t**) group_arrayA, (const int*) group_lda,
		(const BBLAS_Complex64_t**) group_arrayx, (const int*) group_incx,
		(const BBLAS_Complex64_t*) group_beta,
		group_arrayy, (const int*) group_incy,
		(const int) group_count, (const int*) group_size,
		infog);

	printf("Computing with third API...\n");
	// Fixed batch
	enum BBLAS_OPTS batch_opts = BBLAS_FIXED;
	batch_zgemv(
		(const enum BBLAS_TRANS*) &zgemvf_example->trans,
		(const int*) &zgemvf_example->m, (const int*) &zgemvf_example->n,
		(const BBLAS_Complex64_t*) &zgemvf_example->alpha,
		(const BBLAS_Complex64_t**) zgemvf_example->arrayA,
		(const int*) &zgemvf_example->lda,
		(const BBLAS_Complex64_t**) zgemvf_example->arrayx,
		(const int*) &zgemvf_example->incx,
		(const BBLAS_Complex64_t*) &zgemvf_example->beta,
		(BBLAS_Complex64_t**) zgemvf_example->arrayy,
		(const int*) &zgemvf_example->incy,
		(const int) zgemvf_example->batch_count, (const enum BBLAS_OPTS) batch_opts,
		&infof);

	batch_opts = BBLAS_VARIABLE;

	batch_zgemv(
		(const enum BBLAS_TRANS*) zgemvv_example->trans,
		(const int*) zgemvv_example->m, (const int*) zgemvv_example->n,
		(const BBLAS_Complex64_t*) zgemvv_example->alpha,
		(const BBLAS_Complex64_t**) zgemvv_example->arrayA,
		(const int*) zgemvv_example->lda,
		(const BBLAS_Complex64_t **) zgemvv_example->arrayx,
		(const int*) zgemvv_example->incx,
		(const BBLAS_Complex64_t*) zgemvv_example->beta,
		(BBLAS_Complex64_t **) zgemvv_example->arrayy,
		(const int*) zgemvv_example->incy,
		(const int) zgemvv_example->batch_count, (const enum BBLAS_OPTS) batch_opts,
		infov);

	printf("Completed computations with all APIs...\n");

    return 0;
}
