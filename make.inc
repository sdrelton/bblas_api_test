
##############################################################
# BBLAS API TEST make.inc file
# Samuel Relton, Pedro Valero-Lara, Mawussi Zounon
# 5th July 2016
##############################################################

# BBLAS Source Code
BBLAS_BASE_DIR  = $(shell pwd)
BBLAS_UTILS_DIR = $(BBLAS_BASE_DIR)/utils
BBLAS_SRC_DIR   = $(BBLAS_BASE_DIR)/src
BBLAS_TEST_DIR  = $(BBLAS_BASE_DIR)/examples

DEPS            = -I$(BBLAS_BASE_DIR)/include -I$(BBLAS_TEST_DIR)
LDFLAGS         = -fopenmp

# Compiler
CC = gcc
CFLAGS = -c -std=c99 -DADD_ -Wall -pedantic -fopenmp -g
CFLAGS += -DBBLAS_WITH_MKL
DEPS   +=  -m64 -I${MKLROOT}/include
# BLAS libraries
BLAS_LIB =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl

# CBLAS libraries
#CBLAS_DIR       = 
#CBLAS_LIB       = -L$(CBLAS_DIR)/lib
#CBLAS_INC       = -I$(CBLAS_DIR)/include

# LAPACK libraries
#LAPACK_DIR      = 
#LAPACK_LIB      = -L$(LAPACK_DIR)/lib
#LAPACK_INC      = -I$(LAPACK_DIR)/include

# LAPACKE libraries
#LAPACKE_DIR     = 
#LAPACKE_LIB     = -L$(LAPACKE_DIR)/lib -llapacke -llapack
#LAPACKE_INC     = -I$(LAPACKE_DIR)/include

DEPS += $(LAPACKE_INC) $(LAPACK_INC) $(CBLAS_INC)

LDFLAGS += $(LAPACKE_LIB) $(LAPACK_LIB) $(CBLAS_LIB) $(BLAS_LIB) -lm -lgfortran

