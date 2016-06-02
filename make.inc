
##############################################################
# BBLAS make.inc file
#
# Automatically generated by setup.py
#
# BBLAS is a software package provided by Univ. of Manchester,
# Univ. of Tennessee
##############################################################


# BBLAS Source Code
BBLAS_BASE_DIR  = /home/srelton/NLAFET/bblas_api_test
BBLAS_SRC_DIR   = $(BBLAS_BASE_DIR)/src
BBLAS_TEST_DIR  = $(BBLAS_BASE_DIR)/examples

DEPS            = -I$(BBLAS_BASE_DIR)/include -I$(BBLAS_TEST_DIR)
LDFLAGS         =

# Compiler
CC = gcc
CFLAGS = -c -std=c99 -DADD_ -Wall -pedantic
DEPS   += -m64
# BLAS libraries
#BLAS_LIB        = -Wl,--start-group /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64/libmkl_core.a /opt/intel/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl
BLAS_LIB = -lblas

# CBLAS libraries
CBLAS_DIR       = /home/srelton/NLAFET/BBLAS-Manc/install
CBLAS_LIB       = -lcblas # -L$(CBLAS_DIR)/lib -lcblas
CBLAS_INC       = -I$(CBLAS_DIR)/include

# LAPACK libraries
LAPACK_DIR      = /home/srelton/NLAFET/BBLAS-Manc/install
LAPACK_LIB      = -llapack # -L$(LAPACK_DIR)/lib -llapack
LAPACK_INC      = -I$(LAPACK_DIR)/include

# LAPACKE libraries
LAPACKE_DIR     = /home/srelton/NLAFET/BBLAS-Manc/install
LAPACKE_LIB     = -L$(LAPACKE_DIR)/lib -llapacke
LAPACKE_INC     = -I$(LAPACKE_DIR)/include


DEPS += $(LAPACKE_INC) $(LAPACK_INC) $(CBLAS_INC)

LDFLAGS += $(LAPACKE_LIB) $(LAPACK_LIB) $(CBLAS_LIB) $(BLAS_LIB) -ltmglib -lm -lgfortran