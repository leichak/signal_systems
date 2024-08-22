#include "Utils.h"

#include <stdlib.h>

int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

size_t max_size_t(size_t a, size_t b)
{
    return (a > b) ? a : b;
}

void fill_n_with_step(double *xs, size_t size, float start, float stop)
{
    float step = (stop - start) / (float)(size);
    xs[0] = start;
    for (int i = 1; i < size; i++) {
        xs[i] = xs[i - 1] + step;
    }
}
