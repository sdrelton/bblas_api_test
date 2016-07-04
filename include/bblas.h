/**
 * @file bblas.h
 *
 *  @brief BBLAS headers for all precisions.
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
 * Contains all the headers (all precisions) related to batch blas source code
 **/

#ifndef BBLAS_H
#define BBLAS_H

/*
 * BBLAS Headers
 */
#include "bblas_z.h"
//#include "bblas_s.h"
//#include "bblas_d.h"
//#include "bblas_c.h"

/*
 * Error handler
 */

int xerbla_batch(char *func_name, int error, int subproblem_ind);

#endif        //  #ifndef BBLAS_H
