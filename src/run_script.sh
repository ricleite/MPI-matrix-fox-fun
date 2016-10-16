
set -e

TOTAL_PROCS=9

# setup compiler flags
OMPI_CFLAGS="$OMPI_CFLAGS -g -std=c99 -Wall"
export OMPI_CFLAGS

# setup linker flags, must link with math lib
OMPI_LDFLAGS="$OMPI_LDFLAGS -lm"
export OMPI_LDFLAGS

mpicc main.c -o main.out

mpirun -np $TOTAL_PROCS main.out ../tests/input6
