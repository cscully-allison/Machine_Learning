# Produced by CVXGEN, 2018-04-03 18:09:46 -0400.
# CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
# The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
# CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
# applications without prior written permission from Jacob Mattingley.

# Filename: Makefile.
# Description: Basic Makefile.
OPT = #-Wall -Os
# libmath is needed for sqrt, which is used only for reporting the gap. Can
# remove if desired for production solvers..
LDLIBS = -lm
CFLAGS = $(OPT) $(INCLUDES)
CC = nvcc
.PHONY: all
all: test
test: solver.o matrix_support.o ldl.o util.o
# Include util.o for random functions and easy matrix printing.
solver.o: solver.h
matrix_support.o: solver.h
ldl.o: solver.h
util.o: solver.h
.PHONY : clean
clean :
	-rm -f *.o test
