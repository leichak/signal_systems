#include "FixedPoint.h"

/*
 * Convert floating point to fixed point.
 *
 * @param x The floating point number to convert.
 * @param k The number of fractional bits in the fixed point representation.
 * @return The fixed point representation of the input floating point number.
 */
int to_fix(float x, int k)
{
    return (int)(x * (1 << k)); // equivalent to multiplying by 2^k
}

/*
 * Convert fixed point to floating point.
 *
 * @param x The fixed point number to convert.
 * @param k The number of fractional bits in the fixed point representation.
 * @return The floating point representation of the input fixed point number.
 */
float to_float(int x, int k)
{
    return (x * 1.0 / (1 << k)); // equivalent to dividing by 2^k
}

/*
 * Multiply two fixed point numbers.
 *
 * @param x The first fixed point number.
 * @param y The second fixed point number.
 * @param k The number of fractional bits in the fixed point representation.
 * @return The fixed point representation of the product (result is 2N, 2k).
 */
int fixed_N_k_mul(int x, int y, int k)
{
    return (x * y) >> k; // Right shift to adjust the result back to fixed point
}

/*
 * Add two fixed point numbers.
 *
 * @param x The first fixed point number.
 * @param y The second fixed point number.
 * @return The fixed point representation of the sum (result is N+1, k).
 */
int fixed_N_k_add(int x, int y)
{
    return (x + y); // Direct addition, results in potential overflow
}

/*
 * Add two fixed point numbers with different fractional bit representations.
 *
 * @param x_k1 The first fixed point number with k1 fractional bits.
 * @param y_k2 The second fixed point number with k2 fractional bits.
 * @param k1 The number of fractional bits for the first number.
 * @param k2 The number of fractional bits for the second number.
 * @return The fixed point representation of the sum (may result in loss of precision).
 */
int fixed_N1_k2_add(int x_k1, int y_k2, int k1, int k2)
{
    return ((x_k1 >> abs(k1 - k2)) + y_k2); // Right shift x_k1 to align with y_k2
}

/*
 * Add two fixed point numbers with different fractional bit representations.
 *
 * @param x_k1 The first fixed point number with k1 fractional bits.
 * @param y_k2 The second fixed point number with k2 fractional bits.
 * @param k1 The number of fractional bits for the first number.
 * @param k2 The number of fractional bits for the second number.
 * @return The fixed point representation of the sum (may result in overflow).
 */
int fixed_N1_k1_add(int x_k1, int y_k2, int k1, int k2)
{
    return (x_k1 + (y_k2 << (abs(k1 - k2)))); // Left shift y_k2 to align with x_k1
}
