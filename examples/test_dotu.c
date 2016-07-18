/**
 *
 * @file test_dotu.c
 *
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-06-15
 * @precisions normal z -> c
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include<time.h>
#ifdef BBLAS_WITH_MKL
    #include <mkl_cblas.h>
    #include <mkl_lapacke.h>
#else
    #include <cblas.h>
    #include <lapacke.h>
#endif
#include "bblas_z.h"


#undef REAL
#define COMPLEX



int main(int argc, char *argv[])
{
  //Set parameters to initilize random value generation
  __time_t t;
  srand((unsigned) time(&t));
  int IONE     = 1;
  int ISEED[4] ={0,0,0,1};
  
  // Generation of batch_count value between 100 and 1000
  int batch_min = 100;
  int batch_max = 10000;
  int batch_count = rand() % (batch_max - batch_min + 1) + batch_min;

  // Generation of x and y vector sizes
  int * n = (int *) malloc(batch_count*sizeof(int));
  int n_max = 1000;
  int n_min = 10;
  for(int iter=0; iter < batch_count; iter++)
    n[iter] = rand() % (n_max - n_min + 1) + n_min;

  // common calling parameters
  int * incx = (int*) malloc(batch_count*sizeof(int));
  int * incy = (int*) malloc(batch_count*sizeof(int));
  int * info = (int *) malloc(batch_count*sizeof(int));
  enum BBLAS_OPTS batch_opts = BBLAS_VARIABLE;

  //Generation of x and y vectors  and the result vector
  // for using variable api
  BBLAS_Complex64_t ** x, ** y, *result;
  x = (BBLAS_Complex64_t **) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
  y = (BBLAS_Complex64_t **) malloc(batch_count*sizeof(BBLAS_Complex64_t*));
  result = (BBLAS_Complex64_t *) malloc(batch_count*sizeof(BBLAS_Complex64_t));
  
  for(int iter=0; iter < batch_count; iter++){
    x[iter] = (BBLAS_Complex64_t*) malloc(n[iter]*sizeof(BBLAS_Complex64_t));
    y[iter] = (BBLAS_Complex64_t*) malloc(n[iter]*sizeof(BBLAS_Complex64_t));
    
#if defined(BBLAS_WITH_MKL)
    LAPACKE_zlarnv_work(IONE, ISEED, n[iter], (MKL_Complex16*) x[iter]);
    LAPACKE_zlarnv_work(IONE, ISEED, n[iter], (MKL_Complex16*) y[iter]);
#else
    LAPACKE_zlarnv_work(IONE, ISEED, n[iter], x[iter]);
    LAPACKE_zlarnv_work(IONE, ISEED, n[iter], y[iter]);
#endif
    // set values of incx and incy
    incx[iter] = 1;
    incy[iter] = 1;
  }

  //Calling variable api
  printf("Calling variable api: batchv_zdotu_sub \n");
  batchv_zdotu_sub(
		  (const int *) n,
		  (BBLAS_Complex64_t const * const *) x,
		  (const int *) incx,
		  (BBLAS_Complex64_t const * const *) y,
		  (const int *) incy,
		  result,
		  (const int) batch_count, info);
  
  //Calling tag api to solve variable problem
  printf("Calling tag api to solve variable problem: batch_zdotu_sub\n");
  
    batch_zdotu_sub(
		  (const int *) n,
		  (BBLAS_Complex64_t const * const *) x,
		  (const int *) incx,
		  (BBLAS_Complex64_t const * const *) y,
		  (const int *) incy,
		  result,
		  (const int) batch_count,
		  (const enum BBLAS_OPTS) batch_opts,
		  info);

    //Free x and y vectors
    printf("Free x and y vectors after variable case\n");
    
    for(int iter=0; iter < batch_count; iter++){
      free(x[iter]);
      free(y[iter]);
    }
    // free n incx and incy and info
    free(n);
    free(incx);
    free(incy);
    free(info);
    
    // Calling fixe api
    //Generation of x and y vectors, result will be overwritten
    int n_fixe = rand() % (n_max - n_min + 1) + n_min;
    int incx_fixe = 1;
    int incy_fixe = 1;
    int * info_fixe = (int*)malloc(sizeof(int));
    
    for(int iter=0; iter < batch_count; iter++){
      x[iter] = (BBLAS_Complex64_t*) malloc(n_fixe*sizeof(BBLAS_Complex64_t));
      y[iter] = (BBLAS_Complex64_t*) malloc(n_fixe*sizeof(BBLAS_Complex64_t));
      
#if defined(BBLAS_WITH_MKL)
      LAPACKE_zlarnv_work(IONE, ISEED, n_fixe, (MKL_Complex16*) x[iter]);
      LAPACKE_zlarnv_work(IONE, ISEED, n_fixe, (MKL_Complex16*) y[iter]);
#else
      LAPACKE_zlarnv_work(IONE, ISEED, n_fixe, x[iter]);
      LAPACKE_zlarnv_work(IONE, ISEED, n_fixe, y[iter]);
#endif
    }
  //Calling fixe api
    printf("Calling fixe api: batchf_zdotu_sub \n");
  batchf_zdotu_sub(
		  (const int ) n_fixe,
		  (BBLAS_Complex64_t const * const *) x,
		  (const int) incx_fixe,
		  (BBLAS_Complex64_t const * const *) y,
		  (const int) incy_fixe,
		  result,
		  (const int) batch_count, info_fixe);

  //Calling tag api to solve fixed batch  problem
  printf("Calling tag api to solve fixed batch  problem: batch_zdotu_sub \n");
  batch_opts = BBLAS_FIXED;

  batch_zdotu_sub(
		  (const int *) &n_fixe,
		  (BBLAS_Complex64_t const * const *) x,
		  (const int *) &incx_fixe,
		  (BBLAS_Complex64_t const * const *) y,
		  (const int *) &incy_fixe,
		  result,
		  (const int) batch_count,
		  (const enum BBLAS_OPTS) batch_opts,
		  info_fixe);
  //Free x and y vectors
  printf("Free memory after the call to fixed batch case\n");
  
  for(int iter=0; iter < batch_count; iter++){
    free(x[iter]);
    free(y[iter]);
  }

  //Creation of variable specific to group api
  int group_min = 1;
  int group_max = 10;
  int group_count = rand() % (group_max - group_min + 1) + group_min;
  int * group_size = (int *) malloc(group_count*sizeof(int));
  int tmp_size = batch_count/group_count;
  
  //Set value for group_size
  for(int iter=0; iter < group_count-1; iter++)
    group_size[iter] = tmp_size;
  
  //Set last value of group_size
  group_size[group_count-1] = batch_count - (tmp_size*(group_count-1));

  //Allocate memory for vectors

  int *n_group    = (int *) malloc(group_count*sizeof(int));
  int *incx_group = (int *) malloc(group_count*sizeof(int));
  int *incy_group = (int *) malloc(group_count*sizeof(int));
  int *info_group = (int *) malloc(group_count*sizeof(int));
  int offset = 0;
  //Set vector entries

  for(int group_iter = 0; group_iter < group_count ; group_iter++){
    n_group[group_iter] = rand() % (n_max - n_min + 1) + n_min;
    incx_group[group_iter] = 1;
    incy_group[group_iter] = 1;

    for (int iter = 0; iter < group_size[group_iter]; iter++){
      x[offset+iter] = (BBLAS_Complex64_t*) malloc(n_group[group_iter]*sizeof(BBLAS_Complex64_t));
      y[offset+iter] = (BBLAS_Complex64_t*) malloc(n_group[group_iter]*sizeof(BBLAS_Complex64_t));

#if defined(BBLAS_WITH_MKL)
      LAPACKE_zlarnv_work(IONE, ISEED, n_group[group_iter], (MKL_Complex16*) x[offset+iter]);
      LAPACKE_zlarnv_work(IONE, ISEED, n_group[group_iter], (MKL_Complex16*) y[offset+iter]);
#else
      LAPACKE_zlarnv_work(IONE, ISEED, n_group[group_iter], x[offset+iter]);
      LAPACKE_zlarnv_work(IONE, ISEED, n_group[group_iter], y[offset+iter]);
#endif
    }
    offset += group_size[group_iter];
  }

  //Calling group api
  printf("Calling group api:batchg_zdotu_sub\n");
   batchg_zdotu_sub(
		    (const int *)n_group,
		    (BBLAS_Complex64_t const * const *)x,
		    (const int *)incx_group,
		    (BBLAS_Complex64_t const * const *)y,
		    (const int *)incy_group,
		    result,
		    (const int)group_count,
		    (const int *) group_size,
		    info_group);

   //Free allocated memories
  for(int iter=0; iter < batch_count; iter++){
    free(x[iter]);
    free(y[iter]);
  }
  free(result);
  free(x);
  free(y);
  free(n_group);
  free(incx_group);
  free(incy_group);
  free(group_size);
  free(info_group);
  return 0;
}
