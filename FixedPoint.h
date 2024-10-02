#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP

#include <math.h>
#include <stdlib.h>

/*
 * Floating Point to Fixed Point Conversion and Arithmetic
 *
 * This header file covers the conversion between floating point and fixed point
 * representations, including signed and unsigned fixed point formats, scaling,
 * and arithmetic operations. Below are the key concepts and example code for
 * reference.
 *
 * 1. **Increased Complexity of Floating Point**
 *    - **Power Consumption**: Floating point operations can be power-intensive.
 *    - **Energy Efficiency**: Faster programs that use floating point operations
 *      can lead to increased energy consumption.
 *
 * 2. **Fixed Point Refinement**
 *    - **Conversion**: Floating point numbers are converted to fixed-point representation
 *      to optimize power and speed. This process involves scaling and converting floating
 *      point numbers to integer data types.
 *
 * 3. **Unsigned Fixed Point Representation**
 *    - **Representation**: An unsigned integer with bits b_(n-1) to b_0.
 *    - **Value Formula**: V_uint<n> = b_(n-1) * 2^(n-1) + ... + b_0 * 2^0.
 *    - **Fixed Point**: The binary point shifts from the rightmost position to a different
 *      position k. The representation is denoted as ufix<N,k>.
 *    - **Value Formula for ufix<N,k>**:
 *      V_ufix<N,k> = b_(N-k-1) * 2^(N-k-1) + ... + b_0 + b_(-1) * 2^(-1) + ... + b_k * 2^(-k).
 *    - **Example**:
 *
 *      | 2^1 | 2^0 | 2^-1 | 2^-2 | 2^-3 | 2^-4 | Unsigned | ufix<6,4> |
 *      |-----|-----|------|------|------|------|----------|----------|
 *      | 0   | 0   | 0    | 0    | 0    | 0    | 0        | 0        |
 *      | 0   | 0   | 0    | 0    | 0    | 1    | 1        | 0.0625   |
 *      | 0   | 0   | 1    | 0    | 0    | 1    | 9        | 0.5625   |
 *      | 0   | 1   | 0    | 0    | 0    | 0    | 16       | 1.0      |
 *      | 1   | 1   | 1    | 1    | 1    | 1    | 63       | 3.9375   |
 *
 * 4. **Conversion Relationships**
 *    - **Unsigned Integer to Fixed Point**:
 *      - V_uint<N> = V_ufix<N,k> * 2^k
 *      - V_ufix<N,k> = V_uint<N> * 2^(-k)
 *
 *    - **Floating Point to Fixed Point**: Convert by scaling and converting to an
 *      integer data type.
 *
 * 5. **Signed Fixed Point Representation**
 *    - **Two's Complement**: Common method for representing signed integers.
 *    - **Value Formula for Signed Integer**:
 *      V_int<N> = b_(N-1) * [-(2^(N-1))] + ... + b_0.
 *    - **Value Formula for Signed Fixed Point**:
 *      V_fix<N,k> = b_(N-k-1) * [-(2^(N-k-1))] + ... + b_0 + b_(-1) * 2^(-1) + ... + b_k * 2^(-k).
 *    - **Example**:
 *
 *      | 2^1 | 2^0 | 2^-1 | 2^-2 | 2^-3 | 2^-4 | Signed  | ufix<6,4> |
 *      |-----|-----|------|------|------|------|---------|----------|
 *      | 0   | 0   | 0    | 0    | 0    | 0    | 0       | 0        |
 *      | 0   | 0   | 0    | 0    | 0    | 1    | 1       | 0.0625   |
 *      | 1   | 0   | 1    | 0    | 0    | 0    | -24     | -1.5     |
 *      | 1   | 1   | 0    | 0    | 0    | 0    | -16     | -1.0     |
 *      | 1   | 1   | 1    | 1    | 1    | 1    | -1      | -0.0625  |
 *
 * 6. **Conversion Examples**
 *    - **From Floating Point to Fixed Point**:
 *      - Code Example:
 *
 *        ```c
 *        #include <stdio.h>
 *
 *        float taps[5] = {0.1, 1.2, -0.3, 3.4, -1.5};
 *        int taps_8_6[5]; // a fix<8,6> type
 *        unsigned i;
 *
 *        for (i = 0; i < 5; i++)
 *            taps_8_6[i] = (int) (taps[i] * (1 << 6)); // multiply by 2^6
 *
 *        for (i = 0; i < 5; i++)
 *            printf("%+4.2f %10d\n", taps[i], taps_8_6[i]);
 *        ```
 *      - Output:
 *        ```
 *        +0.10          6
 *        +1.20         76
 *        -0.30        -19
 *        +3.40        217
 *        -1.50        -96
 *        ```
 *
 *    - **From Fixed Point to Floating Point**:
 *      - Code Example:
 *
 *        ```c
 *        void fn_fixed_2_float() {
 *            int taps_8_6[5] = {6, 75, -19, 127, -96};
 *            unsigned i;
 *
 *            for (i = 0; i < 5; i++)
 *                taps[i] = (taps_8_6[i] * (1.0 / (1 << 6))); // multiply by 2^-6
 *
 *            for (i = 0; i < 5; i++)
 *                printf("%+4.2f %10d\n", taps[i], taps_8_6[i]);
 *        }
 *        ```
 *
 * 7. **Fixed Point Arithmetic**
 *    - **Unsigned Addition**: Adding two `ufix<N,k>` numbers results in a `ufix<N+1,k>` number,
 *      potentially with overflow or precision loss.
 *    - **Unsigned Subtraction**: Uses the same rules as addition; subtraction can be defined
 *      as the addition of `ufix<N,k>` with the two's complement of another `ufix<N,k>`.
 *    - **Alignment for Addition**: Align `ufix<N,k1>` and `ufix<N,k2>` (k1 > k2) by left or right shifting
 *      before addition, which may result in overflow or precision loss.
 *    - **Multiplication**: Multiplying two `ufix<N,k>` numbers results in `ufix<2N,2k>`. To capture the result
 *      in `ufix<N,k>`, right shift the result by `k` bits, which can also cause precision loss.
 *
 * 8. **Quantization Noise**
 *    - When using fixed point arithmetic, especially in digital filters, precision loss is more
 *      pronounced compared to floating point arithmetic. Errors add to quantization noise, which
 *      appears as broadband noise in useful signals.
 */

int to_fix(float x, int k);
float to_float(int x, int k);
int fixed_N_k_mul(int x, int y, int k);
int fixed_N_k_add(int x, int y);
int fixed_N1_k2_add(int x_k1, int y_k2, int k1, int k2);
int fixed_N1_k1_add(int x_k1, int y_k2, int k1, int k2);

#endif // FIXEDPOINT_HPP
