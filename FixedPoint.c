#include "FixedPoint.h"

int to_fix(float x, int k)
{
    return (int)(x * (1 << k)); // equivalent 2^k
}

float to_float(int x, int k)
{
    return (x * 1.0 / (1 << k)); // equivalent 2^-k
}

int fixed_N_k_mul(int x, int y, int k) // result is 2N, 2K
{
    return (x * y) >> k;
}

int fixed_N_k_add(int x, int y) // N+1 k result
{
    return (x + y);
}

int fixed_N1_k2_add(int x_k1, int y_k2, int k1, int k2) // loss of precision
{                                                       // k1 > k2
    return ((x_k1 >> abs(k1 - k2)) + y_k2);
}

int fixed_N1_k1_add(int x_k1, int y_k2, int k1, int k2) // overflow
{                                                       // k1 > k2
    return (x_k1 + (y_k2 << (abs(k1 - k2))));
}
