#/bin/bash

set -e

TOTAL_PROCS=4

# setup compiler flags
OMPI_CFLAGS="$OMPI_CFLAGS -g -std=c99 -Wall"
export OMPI_CFLAGS

# setup linker flags, must link with math lib
OMPI_LDFLAGS="$OMPI_LDFLAGS -lm"
export OMPI_LDFLAGS

mpicc ../src/main.c -o main.out
for TEST in 6 12 60 300; do
    mpirun -np $TOTAL_PROCS main.out input${TEST} --test output${TEST} --no-matrix-output
done
