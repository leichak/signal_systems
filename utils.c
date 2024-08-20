#include "Utils.h"

int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

void fill_n_with_step(float xs[], float start, float step, int count)
{
    xs[0] = start;
    for (int i = 1; i < count; i++) {
        xs[i] = xs[i - 1] + step;
    }
}
