/**
 * @file xerbla_batch.c
 *
 * @brief BBLAS error handling routine.
 *
 * BBLAS is a software package provided by the Univ. of Tennessee and
 * the Univ. of Manchester
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date 2016-03-04
 *
 **/

#include "bblas_macros.h"
#include <stdio.h>
#include <stdlib.h>

int xerbla_batch(char *func_name, int error, int subproblem_ind)
{
/** Purpose
    -------

    <b>xerbla_batch</b> is an error handler for the BBLAS routines.
    It is called by a BBLAS routine if an input parameter has an
    invalid value. A message is printed.

    Users and vendors may consider modifying this file in order to
    call system-specific exception-handling facilities.


    Parameters
    ----------
	@param[in]
    func_name           <tt>*char</tt>
                        The name of the routine which called xerbla_batch.

	@param[in]
    error               <tt>int</tt>
                        Error code as defined in bblas_macros.h.

	@param[in]
    subproblem_ind      <tt>int</tt>
                        Index of the subproblem where the error orginated.
**/

    switch(error)
    {
        case BBLAS_ERR_BATCH_COUNT:
            fprintf(stderr, "BBLAS ERROR in function %s caused by BATCH_COUNT.\n",func_name);
	    return(EXIT_SUCCESS);
            /* To terminate replace "break;" with "exit(EXIT_FAILURE);" */
            break;
        case BBLAS_ERR_BATCH_OPTS:
            fprintf(stderr, "BBLAS ERROR in function %s caused by BATCH_OPTS.\n",func_name);
	    return(EXIT_SUCCESS);
		    /* To terminate replace "break;" with "exit(EXIT_FAILURE);" */
            break;
        case BBLAS_ERR_M:
            fprintf(stderr, "BBLAS ERROR in function %s caused by M[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_N:
            fprintf(stderr, "BBLAS ERROR in function %s caused by N[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_K:
            fprintf(stderr, "BBLAS ERROR in function %s caused by K[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_LDA:
            fprintf(stderr, "BBLAS ERROR in function %s caused by LDA[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_LDB:
            fprintf(stderr, "BBLAS ERROR in function %s caused by LDB[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_LDC:
            fprintf(stderr, "BBLAS ERROR in function %s caused by LDC[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_UPLO:
            fprintf(stderr, "BBLAS ERROR in function %s caused by UPLO[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_TRANSA:
            fprintf(stderr, "BBLAS ERROR in function %s caused by TRANSA[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_TRANSB:
            fprintf(stderr, "BBLAS ERROR in function %s caused by TRANSB[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_TRANS:
            fprintf(stderr, "BBLAS ERROR in function %s caused by TRANS[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_SIDE:
            fprintf(stderr, "BBLAS ERROR in function %s caused by SIDE[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        case BBLAS_ERR_DIAG:
            fprintf(stderr, "BBLAS ERROR in function %s caused by DIAG[%d].\n",func_name,subproblem_ind);
            return(0);
            break;
        default:
            fprintf(stderr, "BBLAS ERROR in function %s caused by undefined error code.\n",func_name);
            exit(EXIT_FAILURE);
            break;
    }
} /* End of xerbla_batch */
