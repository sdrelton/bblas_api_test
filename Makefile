include ./make.inc

BBLAS_SRC_LIST   = batchf_zgemv.c batchv_zgemv.c batchg_zgemv.c batch_zgemv.c \
				   batchf_zgemm.c batchv_zgemm.c batchg_zgemm.c batch_zgemm.c \
		           xerbla_batch.c

BBLAS_SRC=$(addprefix $(BBLAS_SRC_DIR)/, $(BBLAS_SRC_LIST))

TEST_SRC      =  random_gen.c
TEST_SRC_LIST=$(addprefix $(BBLAS_TEST_DIR)/, $(TEST_SRC))

SOURCES       = $(BBLAS_SRC) $(TEST_SRC_LIST)
SOURCES_Z      = $(SOURCES)
OBJECTS_Z     = $(SOURCES_Z:.c=.o)

all:
	$(MAKE) testzgemv
	$(MAKE) testzgemm

.DEFAULT_GOAL := all

testzgemv: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_gemv.c -o $(BBLAS_TEST_DIR)/test_gemv.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_gemv.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

testzgemm: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_gemm.c -o $(BBLAS_TEST_DIR)/test_gemm.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_gemm.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

.c.o:
	$(CC) $(CFLAGS) $(DEPS) $<   -o $@

clean:
	rm */*.o
	rm */testzgemv
	rm */testzgemm
