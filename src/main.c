
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

#define MATRIX_INFINITY 200000

void print_usage();
int main_impl(int argc, char** argv, int num_tasks, int world_rank);
int read_matrix(char* file_str, Matrix* m);
void compute_matrix(Matrix* master_m, Matrix* m, int world_rank,
    MPI_Comm grid_comm, MPI_Comm row_comm, MPI_Comm col_comm);
void scatter_matrix_to_slaves(Matrix const* master_m, Matrix* m, int world_rank);
void gather_matrix_from_slaves(Matrix* master_m, Matrix* m, int world_rank);
void convert_matrix(Matrix* m, int for_output);
void print_matrix(Matrix const* m);

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int num_tasks;
    int world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // MPI_Finalize has to be always called, so wrap error checks/functionality
    int exit_code = main_impl(argc, argv, num_tasks, world_rank);
    if (exit_code > 1 && world_rank == RANK_MASTER)
        print_usage();

    MPI_Finalize();
    return exit_code;
}

void print_usage()
{
    printf("Usage: ./main.out matrix_file [options]\n");
    printf("Options:\n");
    printf("'--test test_file' uses matrix in test_file to compare with output matrix, exit error 1 if not equal\n");
    printf("'--no-matrix-output' silences default matrix print output\n");
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

    char* test_file_str = NULL;
    int matrix_output = 2;

    // master has to:
    // - process program options
    // - do error checking, notify slaves to exit early if error
    // - read the entire matrix file and distribute it
    if (world_rank == RANK_MASTER)
    {
        if (argc <= 1) // missing matrix file arg
            master_exit_code = 2;
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
                    master_exit_code = 2;
                }
                else
                {
                    // now we know slice size
                    m.n = master_m.n / q;
                }
            }

            // process program options
            if (master_exit_code == 0)
            {
                for (int i = 2; i < argc; ++i)
                {
                    char* arg = argv[i];
                    if (strcmp(arg, "--test") == 0)
                    {
                        ++i;
                        if (i >= argc)
                        {
                            LOG_ERROR("Missing test_file arg");
                            master_exit_code = 2;
                            break;
                        }

                        test_file_str = argv[i];
                    }
                    else if (strcmp(arg, "--no-matrix-output") == 0)
                        matrix_output = 0;
                    else
                    {
                        LOG_ERROR("Invalid arg '%s'", arg);
                        master_exit_code = 2;
                        break;
                    }
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

    // setup grid communicators
    MPI_Comm grid_comm;
    MPI_Comm row_comm;
    MPI_Comm col_comm;

    {
        int slices_per_line = master_m.n / m.n;
        int dims[2] = { slices_per_line, slices_per_line };
        int periods[2] = { 1, 1 };
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &grid_comm);

        int row_dims[2] = { 0, 1 };
        MPI_Cart_sub(grid_comm, row_dims, &row_comm);

        int col_dims[2] = { 1, 0 };
        MPI_Cart_sub(grid_comm, col_dims, &col_comm);
    }

    compute_matrix(&master_m, &m, world_rank, grid_comm, row_comm, col_comm);

    if (world_rank == RANK_MASTER)
    {
        if (matrix_output)
            print_matrix(&master_m);

        if (test_file_str)
        {
            Matrix test;
            Matrix_Init(&test);
            if (read_matrix(test_file_str, &test))
                master_exit_code = 2; // failed to read matrix
            else
            {
                master_exit_code = Matrix_Compare(&master_m, &test);
                if (master_exit_code == 0)
                    printf("Passed test for matrix size %d\n", master_m.n);
                else
                    printf("Failed test for matrix size %d\n", master_m.n);
            }
        }
    }

    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&grid_comm);

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

void compute_matrix(Matrix* master_m, Matrix* m, int world_rank,
    MPI_Comm grid_comm, MPI_Comm row_comm, MPI_Comm col_comm)
{
    // replace 0s by infinite cost, 0 cost between 2 nodes means no link
    if (world_rank == RANK_MASTER)
        convert_matrix(master_m, 0);

    // perform Floyd-Warshall to compute shortest paths
    int const iterations = (int)ceil(log2(master_m->n));
    int const slice_size = m->n;
    int const slices_per_line = master_m->n / slice_size;

    assert(pow(2, iterations) > master_m->n);

    int grid_rank;
    MPI_Comm_rank(grid_comm, &grid_rank);

    int coords[2];
    MPI_Cart_coords(grid_comm, grid_rank, 2, coords);

    int my_row = coords[0];
    int my_col = coords[1];

    Matrix temp_a;
    Matrix temp_b;
    Matrix temp_c;
    Matrix_Init(&temp_a);
    Matrix_Init(&temp_b);
    Matrix_Init(&temp_c);
    Matrix_Create(&temp_a, slice_size);
    Matrix_Create(&temp_b, slice_size);
    Matrix_Create(&temp_c, slice_size);

    // need to run log(n) matrix multiplications for Floyd-Warshall algorithm
    for (int itr = 0; itr < iterations; ++itr)
    {
        scatter_matrix_to_slaves(master_m, m, world_rank);

        // run Fox's algorithm
        // A = temp_a, B = temp_b, C = temp_c, original matrix = m
        memcpy(temp_b.mat, m->mat, (slice_size * slice_size) * sizeof(int));
        memset(temp_c.mat, MATRIX_INFINITY, (slice_size * slice_size) * sizeof(int));

        for (int stage = 0; stage < slices_per_line; ++stage)
        {
            int cast_root = (my_row + stage) % slices_per_line;
            Matrix* a = (cast_root == my_col) ? m : &temp_a;
            Matrix* b = &temp_b;
            Matrix* c = &temp_c;

            MPI_Bcast(a->mat, slice_size * slice_size, MPI_INT, cast_root, row_comm);

            for (int i = 0; i < slice_size; ++i)
            {
                for (int j = 0; j < slice_size; ++j)
                {
                    int c_value = Matrix_GetValue(c, i, j);
                    for (int row_col = 0; row_col < slice_size; ++row_col)
                    {
                        int a_value = Matrix_GetValue(a, i, row_col);
                        int b_value = Matrix_GetValue(b, row_col, j);

                        c_value = MIN(c_value, a_value + b_value);
                    }

                    Matrix_SetValue(c, i, j, c_value);
                }
            }

            // dest for our current B
            int const above_in_col = (slices_per_line + my_row - 1) % slices_per_line;
            // source of our next B
            int const below_in_col = (my_row + 1) % slices_per_line;

            MPI_Sendrecv_replace(b->mat, slice_size * slice_size, MPI_INT,
                above_in_col, 0, below_in_col, 0, col_comm, NULL);
        }

        gather_matrix_from_slaves(master_m, &temp_c, world_rank);
    }

    // convert infinite costs back to 0, for output/comparisons
    if (world_rank == RANK_MASTER)
        convert_matrix(master_m, 1);

    Matrix_Delete(&temp_a);
    Matrix_Delete(&temp_b);
    Matrix_Delete(&temp_c);
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
    MPI_Scatter(slice_ptr, (slice_n * slice_n), MPI_INT,
        slice_recv_ptr, (slice_n * slice_n), MPI_INT, RANK_MASTER, MPI_COMM_WORLD);

    if (world_rank == RANK_MASTER)
        free(slice_ptr);
}

void gather_matrix_from_slaves(Matrix* master_m, Matrix* m, int world_rank)
{
    int const matrix_n = master_m->n;
    int const slice_n = m->n;

    int* slice_ptr = NULL;
    int* slice_send_ptr = m->mat;
    if (world_rank == RANK_MASTER)
        slice_ptr = (int*)malloc(matrix_n * matrix_n * sizeof(int));

    MPI_Gather(slice_send_ptr, (slice_n * slice_n), MPI_INT,
        slice_ptr, (slice_n * slice_n), MPI_INT, RANK_MASTER, MPI_COMM_WORLD);

    if (world_rank == RANK_MASTER)
    {
        Matrix_Join(master_m, slice_ptr, slice_n);
        free(slice_ptr);
    }
}

void convert_matrix(Matrix* m, int for_output)
{
    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            if (i == j)
                continue;

            int value = Matrix_GetValue(m, i, j);
            if (for_output && value >= MATRIX_INFINITY)
                Matrix_SetValue(m, i, j, 0);
            else if (for_output == 0 && value == 0)
                Matrix_SetValue(m, i, j, MATRIX_INFINITY);
        }
    }


}

void print_matrix(Matrix const* m)
{
    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            int value = Matrix_GetValue(m, i, j);
            printf("%d ", value);
        }

        printf("\n");
    }
}
