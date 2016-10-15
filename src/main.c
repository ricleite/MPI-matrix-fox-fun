
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "mpi.h"
#include "defines.h"
#include "matrix.h"

#if __STDC_VERSION__ < 199901L
#error "Compile with at least C99 (-std=c99)"
#endif

void print_usage();
int main_impl(int argc, char** argv);
int read_matrix(char* file_str, Matrix* m);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    // MPI_Finalize has to be always called, so wrap error checks/functionality
    // only master returns error codes, to end process early
    int error_code = main_impl(argc, argv);
    if (error_code > 0)
    {
        print_usage();
        MPI_Abort(MPI_COMM_WORLD, error_code);
    }

    MPI_Finalize();
    return error_code;
}

void print_usage()
{
    printf("Usage: ./main.out matrix_file [options]\n");
    printf("Options:\n");
    printf("'--test test_file' uses matrix in test_file to compare with output matrix, exit error 1 if not equal");
    printf("'--help' prints this usage message");
}

int main_impl(int argc, char** argv)
{
    int num_tasks;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // master has to:
    // - do error checking, notify slaves to exit immediately if error
    // - read the entire matrix file and distribute it
    if (rank == RANK_MASTER)
    {
        if (argc <= 1)
            return 1;

        char* in_file_str = argv[1];

        Matrix m;
        int error = read_matrix(in_file_str, &m);
        if (error)
            return 1;

        // @todo
    }
    else
    {
        // @todo
    }

    printf("num_tasks %d, rank %d\n", num_tasks, rank);
    return 0;
}

int read_matrix(char* file_str, Matrix* m)
{
    assert(file_str);
    assert(m);

    FILE* file = fopen(file_str, "r");
    if (!file)
    {
        LOG_ERROR("couldn't open file %s", file);
        return 1;
    }

    int n;
    if (fscanf(file, "%d ", &n) == 0)
    {
        LOG_ERROR("couldn't read matrix size from file %s", file);
        return 2;
    }

    Matrix_Init(m, n);

    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            int value;
            if (fscanf(file, "%d ", &value) == 0)
            {
                LOG_ERROR("couldn't read matrix value in position (%d, %d), file %s",
                    i, j, file);

                return 4;
            }

            Matrix_SetValue(m, i, j, value);
        }
    }

    return 0;
}
