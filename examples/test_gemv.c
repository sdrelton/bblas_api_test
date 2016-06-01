#include "bblas_z.h"
#include "example_z.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    printf("Testing gemv using some very simple tests.\n");

	/* Create simple test matrices */
	printf("Creating inputs for testing grouped zgemv...\n");

	struct zgemv_batchf_example *zgemv_example =
		(struct zgemv_batchf_example*) malloc(sizeof(struct zgemv_batchf_example));
	set_params_fixed_zgemv(zgemv_example);

	int *info = malloc(sizeof(int)*zgemv_example->batch_count);

	batchf_zgemv(
		(const enum BBLAS_TRANS) zgemv_example->trans,
		(const int) zgemv_example->m, (const int) zgemv_example->n,
		(const BBLAS_Complex64_t) zgemv_example->alpha,
		(const BBLAS_Complex64_t**) zgemv_example->arrayA,
		(const int) zgemv_example->lda,
		(const BBLAS_Complex64_t **) zgemv_example->arrayx,
		(const int) zgemv_example->incx,
		(const BBLAS_Complex64_t) zgemv_example->beta,
		(BBLAS_Complex64_t **) zgemv_example->arrayy,
		(const int) zgemv_example->incy,
		(const int) zgemv_example->batch_count, info);

	printf("Completed fixed computation...\n");

    return 0;
}
