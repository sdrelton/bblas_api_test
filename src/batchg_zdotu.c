/**
 * @file batchg_zdotu.c
 *
 * Part of API test for Batched BLAS routines.
 *
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date 
 *
 * @precisions normal z -> c d s
 *
 **/

#include <cblas.h>
#include "bblas.h"

#define COMPLEX

void batchg_zdotu_sub(
		  const int *n,
		  BBLAS_Complex64_t const * const * x,
		  const int *incx,
		  BBLAS_Complex64_t const * const * y,
		  const int *incy,
		  BBLAS_Complex64_t *dotu,
		  const int group_count, const int *group_size,
		  int* info)
{
	/* Local variables */
	char func_name[15] = "batchg_zdotu";
	int group_iter = 0;
	int offset = 0; // How many subproblems solved so far
	int end_of_group = 0; // End of the current group

	/* Check group_count */
	if (group_count < 0)
	  {
	    xerbla_batch(func_name, BBLAS_ERR_GROUP_COUNT, -1);
	    return;
	  }

	/* Check group_size and call fixed batch computation */
	for (group_iter = 0; group_iter < group_count; group_iter++)
	{
		if (group_size[group_iter] < 0)
		{
			end_of_group = end_of_group + group_size[group_iter];
			xerbla_batch(func_name, BBLAS_ERR_GROUP_SIZE, group_iter);
			info[group_iter] = BBLAS_ERR_GROUP_SIZE;
			continue;
		}

		/* Call fixed batch computation on the group */
		batchf_zdotu_sub(
				 n[group_iter],
				 x+offset,
				 incx[group_iter],
				 y+offset,
				 incy[group_iter],
				 dotu+offset,
				 group_size[group_iter],
				 info+group_iter);
		/* Update offset */
		offset += group_size[group_iter];
	} // End of group loop
}
#undef COMPLEX
