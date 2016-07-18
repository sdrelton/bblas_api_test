include ./make.inc

BBLAS_UTILS      = $(BBLAS_UTILS_DIR)/xerbla_batch.c
BBLAS_SRC_LIST   = batchf_zgemv.c batchv_zgemv.c batchg_zgemv.c batch_zgemv.c \
				   batchf_zgemm.c batchv_zgemm.c batchg_zgemm.c batch_zgemm.c batchf_zgemm_stride.c \
			       batchf_zdotu.c batchv_zdotu.c batchg_zdotu.c batch_zdotu.c

BBLAS_SRC=$(addprefix $(BBLAS_SRC_DIR)/, $(BBLAS_SRC_LIST))

TEST_SRC      = testing_utils.c
TEST_SRC_LIST = $(addprefix $(BBLAS_TEST_DIR)/, $(TEST_SRC))

SOURCES       = $(BBLAS_UTILS) $(BBLAS_SRC) $(TEST_SRC_LIST)
SOURCES_Z     = $(SOURCES)
OBJECTS_Z     = $(SOURCES_Z:.c=.o)

all:
	$(MAKE) test_gemv
	$(MAKE) test_gemm
	$(MAKE) test_dotu
.DEFAULT_GOAL := all

test_dotu: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_dotu.c -o $(BBLAS_TEST_DIR)/test_dotu.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_dotu.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

test_gemv: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_gemv.c -o $(BBLAS_TEST_DIR)/test_gemv.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_gemv.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

test_gemm: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_gemm.c -o $(BBLAS_TEST_DIR)/test_gemm.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_gemm.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

.c.o:
	$(CC) $(CFLAGS) $(DEPS) $<   -o $@

clean:
	rm */*.o
	rm */test_gemv
	rm */test_gemm
	rm */test_dotu
