# Point to includes of local directories
CC=/usr/local/cuda-9.0/bin/nvcc
INCLUDES=-I../include
O_FILES=test.o matrix_support.o solver.o ldl.o

all: $(O_FILES)
	$(CC) -o BuildTest $(O_FILES)

test.o: ../src/test.cu
	$(CC) -dc ../src/test.cu -o test.o $(INCLUDES)

#util.o: ../src/util.cu
#	$(CC) -dc ../src/util.cu -o util.o $(INCLUDES)

matrix_support.o: ../src/matrix_support.cu
	$(CC) -dc ../src/matrix_support.cu -o matrix_support.o $(INCLUDES)

solver.o: ../src/solver.cu
	$(CC) -dc ../src/solver.cu -o solver.o $(INCLUDES)

ldl.o: ../src/ldl.cu
	$(CC) -dc ../src/ldl.cu -o ldl.o $(INCLUDES)

clean:
	rm -f *.o
	rm -f *~
