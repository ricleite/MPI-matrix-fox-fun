
#include <assert.h>
#include <string.h>
#include <unistd.h>

struct Matrix
{
    int n;
    int* mat; // size n*n
};

typedef struct Matrix Matrix;

void Matrix_Init(Matrix* m);
void Matrix_Create(Matrix* m, int n);
void Matrix_Delete(Matrix* m);
void Matrix_Zero(Matrix* m);
int Matrix_GetValue(Matrix const* m, int i, int j);
void Matrix_SetValue(Matrix* m, int i, int j, int value);
int Matrix_Compare(Matrix const* m1, Matrix const* m2);
void Matrix_Split(Matrix const* m, int* slice_ptr, int slice_size);
void Matrix_Join(Matrix* m, int* slice_ptr, int slice_size);

inline void Matrix_Init(Matrix* m)
{
    m->n = 0;
    m->mat = NULL;
}

inline void Matrix_Create(Matrix* m, int n)
{
    assert(m);
    assert(n > 0);

    m->n = n;
    m->mat = (int*)malloc(n * n * sizeof(int));
    Matrix_Zero(m);
}

inline void Matrix_Delete(Matrix* m)
{
    if (m->mat)
        free(m->mat);

    m->mat = NULL;
}

inline void Matrix_Zero(Matrix* m)
{
    memset(m->mat, 0, m->n * sizeof(int));
}

inline int Matrix_GetValue(Matrix const* m, int i, int j)
{
    assert(m);
    assert(m->mat);
    assert(i < m->n && j < m->n);
    return m->mat[(m->n * i) + j];
}

inline void Matrix_SetValue(Matrix* m, int i, int j, int value)
{
    assert(m);
    assert(m->mat);
    assert(i < m->n && j < m->n);
    m->mat[(m->n * i) + j] = value;
}

inline int Matrix_Compare(Matrix const* m1, Matrix const* m2)
{
    assert(m1 && m2);
    assert(m1->mat && m2->mat);
    if (m1->n != m2->n)
        return 1;

    for (int i = 0; i < m1->n; ++i)
    {
        for (int j = 0; j < m2->n; ++j)
        {
            int v1 = Matrix_GetValue(m1, i, j);
            int v2 = Matrix_GetValue(m2, i, j);
            if (v1 != v2)
                return 1;
        }
    }

    return 0;
}

// split matrix[n][n] into ordered slice[slice_size][slice_size] pieces
void Matrix_Split(Matrix const* m, int* slice_ptr, int slice_size)
{
    assert(m);
    assert(m->mat);
    assert(slice_ptr);
    assert(m->n % slice_size == 0);

    int const slices_per_line = m->n / slice_size;

    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            // slice number in the matrix
            int slice_num = (i / slice_size) * slices_per_line + (j / slice_size);
            // calc coords inside slice
            int slice_i = i % slice_size;
            int slice_j = j % slice_size;

            int idx = slice_num * (slice_size * slice_size);
            // add local slice coords
            idx += (slice_i * slice_size) + slice_j;

            int value = Matrix_GetValue(m, i, j);
            slice_ptr[idx] = value;
        }
    }
}

// join ordered slice[slice_size][slice_size] pieces back into a matrix[n][n]
void Matrix_Join(Matrix* m, int* slice_ptr, int slice_size)
{
    assert(m);
    assert(m->mat);
    assert(slice_ptr);
    assert(m->n % slice_size == 0);

    int const slices_per_line = m->n / slice_size;

    for (int i = 0; i < m->n; ++i)
    {
        for (int j = 0; j < m->n; ++j)
        {
            // slice number in the matrix
            int slice_num = (i / slice_size) * slices_per_line + (j / slice_size);
            // calc coords inside slice
            int slice_i = i % slice_size;
            int slice_j = j % slice_size;

            int idx = slice_num * (slice_size * slice_size);
            // add local slice coords
            idx += (slice_i * slice_size) + slice_j;

            int value = slice_ptr[idx];
            Matrix_SetValue(m, i, j, value);
        }
    }
}
