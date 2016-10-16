
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "defines.h"
#include "matrix.h"

#if __STDC_VERSION__ < 199901L
#error "Compile with at least C99 (-std=c99)"
#endif

void print_usage();
int main_impl(int argc, char** argv, int num_tasks, int world_rank);
int read_matrix(char* file_str, Matrix* m);
void scatter_matrix_to_slaves(Matrix const* master_m, Matrix* m, int world_rank);
void gather_matrix_from_slaves(Matrix* master_m, Matrix* m, int world_rank);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int num_tasks;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // MPI_Finalize has to be always called, so wrap error checks/functionality
    int exit_code = main_impl(argc, argv, num_tasks, world_rank);
    if (exit_code && world_rank == RANK_MASTER)
        print_usage();

    MPI_Finalize();
    return exit_code;
}

void print_usage()
{
    printf("Usage: ./main.out matrix_file [options]\n");
    printf("Options:\n");
    printf("'--test test_file' uses matrix in test_file to compare with output matrix, exit error 1 if not equal\n");
    printf("'--help' prints this usage message\n");
}

int main_impl(int argc, char** argv, int num_tasks, int world_rank)
{
    int master_exit_code = 0;

    // full matrix is kept by master, but Matrix::n is used in all slaves
    Matrix master_m;
    Matrix m; // slave matrix slice

    Matrix_Init(&master_m);
    Matrix_Init(&m);

    // master has to:
    // - do error checking, notify slaves to exit early if error
    // - read the entire matrix file and distribute it
    if (world_rank == RANK_MASTER)
    {
        if (argc <= 1) // missing matrix file arg
            master_exit_code = 1;
        else
        {
            char* in_file_str = argv[1];

            if (read_matrix(in_file_str, &master_m))
                master_exit_code = 1; // failed to read matrix
            else
            {
                int q = (int)sqrt(num_tasks);
                // number of slaves needs to be a perfect square
                // that perfect square must divide the matrix size
                if (q * q != num_tasks || master_m.n % q != 0)
                {
                    LOG_ERROR("Bad slave count (%d), can't split matrix (size %d)",
                        num_tasks, master_m.n);
                    master_exit_code = 1;
                }
                else
                {
                    // now we know slice size
                    m.n = master_m.n / q;
                }
            }
        }
    }


    // exit early if master has something to say after error checking
    MPI_Bcast(&master_exit_code, 1, MPI_INT, RANK_MASTER, MPI_COMM_WORLD);
    if (master_exit_code)
        goto main_impl_cleanup;

    // @todo: optimize into a single send
    MPI_Bcast(&master_m.n, 1, MPI_INT, RANK_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&m.n, 1, MPI_INT, RANK_MASTER, MPI_COMM_WORLD);

    // create slave slice for all processes
    Matrix_Create(&m, m.n);

    scatter_matrix_to_slaves(&master_m, &m, world_rank);

    printf("num_tasks %d, world_rank %d\n", num_tasks, world_rank);

main_impl_cleanup:
    if (world_rank == RANK_MASTER)
        Matrix_Delete(&master_m);

    Matrix_Delete(&m);
    return (world_rank == RANK_MASTER) ? master_exit_code : 0;
}

int read_matrix(char* file_str, Matrix* m)
{
    assert(file_str);
    assert(m);

    FILE* file = fopen(file_str, "r");
    if (!file)
    {
        LOG_ERROR("couldn't open file %s", file_str);
        return 1;
    }

    int n;
    if (fscanf(file, "%d ", &n) == 0)
    {
        LOG_ERROR("couldn't read matrix size from file %s", file_str);
        return 2;
    }

    Matrix_Create(m, n);

    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            int value;
            if (fscanf(file, "%d ", &value) == 0)
            {
                LOG_ERROR("couldn't read matrix value in position (%d, %d), file %s",
                    i, j, file_str);

                return 4;
            }

            Matrix_SetValue(m, i, j, value);
        }
    }

    return 0;
}

void scatter_matrix_to_slaves(Matrix const* master_m, Matrix* m, int world_rank)
{
    int const matrix_n = master_m->n;
    int const slice_n = m->n;

    int* slice_ptr = NULL;
    int* slice_recv_ptr = m->mat;
    if (world_rank == RANK_MASTER)
    {
        slice_ptr = (int*)malloc(matrix_n * matrix_n * sizeof(int));
        Matrix_Split(master_m, slice_ptr, slice_n);
    }

    // scatter matrix slices, each one with size (slice_n * slice_n)
    MPI_Scatter((void*)slice_ptr, (slice_n * slice_n), MPI_INT,
        (void*)slice_recv_ptr, (slice_n * slice_n), MPI_INT, RANK_MASTER, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == RANK_MASTER)
        free(slice_ptr);
}

void gather_matrix_from_slaves(Matrix* master_m, Matrix* m, int world_rank)
{

}
