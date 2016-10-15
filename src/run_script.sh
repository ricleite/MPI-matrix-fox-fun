
set -e

# setup compiler flags
OMPI_CFLAGS="$OMPI_CFLAGS -g -std=c99"
export OMPI_CFLAGS

mpicc main.c -o main.out

mpirun -np 5 main.out ../tests/input300
