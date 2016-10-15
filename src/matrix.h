
#include <assert.h>

struct Matrix
{
    int n;
    int* mat; // size n*n
};

typedef struct Matrix Matrix;

void Matrix_Init(Matrix* m, int n);
void Matrix_Delete(Matrix* m);
int Matrix_GetValue(Matrix* m, int i, int j);
void Matrix_SetValue(Matrix* m, int i, int j, int value);
int Matrix_Compare(Matrix* m1, Matrix* m2);

inline void Matrix_Init(Matrix* m, int n)
{
    assert(m);
    assert(n > 0);

    m->n = n;
    m->mat = (int*)malloc(n * n * sizeof(int));
}

inline void Matrix_Delete(Matrix* m)
{
    free(m->mat);
    m->mat = NULL;
}

inline int Matrix_GetValue(Matrix* m, int i, int j)
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

inline int Matrix_Compare(Matrix* m1, Matrix* m2)
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
