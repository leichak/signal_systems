#ifndef FILTERING_H
#define FILTERING_H

#include "DigitalFilters.h"
#include "FixedPoint.h"

/*
 * The Four Direct Forms (https://www.dsprelated.com/freebooks/filters/Four_Direct_Forms.html)
 * =============================================================================
 *
 * Direct-form I (DF1)
 * -------------------
 * Properties:
 *  - DF1 can be viewed as a two-zero filter section followed by a two-pole filter section in series.
 *  - In most fixed-point arithmetic schemes (e.g., two's complement), there is no risk of internal filter overflow
 *    because DF1 has only one summation point. Fixed-point overflow naturally wraps around from the largest positive
 *    value to the largest negative value (and vice versa). As long as the final output y(n) is within range, overflow
 *    is avoided, even if intermediate summations overflow.
 *  - DF1 requires twice as many delays as necessary and is not canonical with respect to delay. However, it is ALWAYS
 *    possible to implement an nth-order filter with N delays.
 *  - The poles and zeros of the filter can be very sensitive to round-off errors in the coefficients. While this isn't
 *    a major issue for simple second-order sections, it can be problematic for higher-order sections.
 *  - The sensitivity to round-off errors is similar to the sensitivity of polynomial roots to coefficient rounding.
 *    This sensitivity tends to be more significant when roots are closely clustered in the complex plane, as opposed
 *    to when they are well spread.
 *  - To minimize this sensitivity, it's common practice to factor the filter transfer function into series or parallel
 *    second-order sections.
 *
 * Two-Complement Wrap-Around
 * --------------------------
 * Temporary overflow in two's complement causes no ill effects.
 * In 3-bit signed fixed-point arithmetic, the available numbers are:
 *
 *     -4  100
 *     -3  101
 *     -2  110
 *     -1  111
 *      0  000
 *      1  001
 *      2  010
 *      3  011
 *
 * Example:
 *  Sum of 3 + 3 - 4 = 2, which gives temporary overflow (3 + 3 = 6 -> wraps to -2),
 *  but the final result (2) is within the allowed range [-4, 3].
 *
 *     011 + 011 = 110  (3 + 3 -> -2)
 *     110 + 100 = 010  (-2 + (-4) -> 2)
 *
 * Similarly, 1 + 3 - 2 = 2:
 *
 *     001 + 011 = 100  (1 + 3 -> -4)
 *     100 + 110 = 010  (-4 + (-2) -> 2)
 *
 * Direct-form II (DF2)
 * -------------------
 * The difference equation for a second-order DF-II structure can be written as:
 *
 *     v(n) = x(n) - a1 * v(n-1) - a2 * v(n-2)
 *     y(n) = b0 * v(n) + b1 * v(n-1) + b2 * v(n-2)
 *
 * This can be interpreted as a two-pole filter followed by a two-zero filter in series. This contrasts with the DF1
 * structure, where the two-zero FIR section precedes the two-pole recursive section in series.
 * Since LTI filters in series commute, we can reverse this ordering and implement an all-pole filter followed by an
 * FIR filter in series. In this configuration, delay elements in both sections contain the same values, allowing
 * a single delay line to be shared between the all-pole and all-zero FIR sections.
 *
 * Properties:
 *  - Can be regarded as a two-pole filter section followed by a two-zero filter section.
 *  - Canonical with respect to delay, as delays are shared.
 *  - In fixed-point arithmetic, overflow can occur at the delay-line input, unlike in DF1.
 *  - Poles and zeros are sensitive to round-off errors in the coefficients a1 and b1, especially in higher-order filters.
 *
 * In DF2, the signal entering the delay line typically requires a larger dynamic range than the output signal y(n).
 * The feedback section of DF2 often amplifies the signal, which is then attenuated by the feedforward (zero) section.
 * To prevent overflow, the two delay elements may need extra guard bits to accommodate the extended dynamic range.
 * Even doubling the number of bits in the delay elements (which doesn't guarantee elimination of internal overflow)
 * can cancel out the memory savings of DF2 over DF1. Thus, while DF2 is canonical with respect to delay, it may require
 * as much or more memory as DF1, even though DF1 uses twice as many delay elements for filter state memory.
 *
 * Transposed Direct Forms
 * -----------------------
 * The remaining two direct forms are obtained by transposing DF1 and DF2. This process is sometimes called
 * "graph reversal." Transposing a single-input, single-output (SISO) filter does not alter its transfer function.
 * This fact can be derived from Mason's gain formula for signal flow graphs or from Tellegen's theorem in network theory.
 *
 * The transpose of a SISO digital filter is straightforward:
 *  - Reverse the direction of all signal paths.
 *  - Make necessary adjustments:
 *      - Convert signal branch-points to summing junctions.
 *      - Convert summing junctions to branch points.
 *      - This will place the input on the right and the output on the left.
 *      - The diagram is then typically flipped horizontally.
 *
 * Numerical Robustness of TDF-II
 * ------------------------------
 * One advantage of the DF2 structure is that zeros precede poles in series order. In many digital filter designs, the
 * poles produce large gains at certain frequencies, which are often compensated for by the zeros. This ordering can
 * improve numerical robustness in some cases.
 */

/// Convention - ak is poles, bk vector is zeros (numerator)
/// coefficients are in order a0, a-1, a-2 ...

/// @brief  Direct-form I
typedef struct
{
    double *b_k;
    double *a_k;
    double *y_k;
    double *x_k;
    size_t size_b;
    size_t size_a;
    size_t max_size;
} DirectForm1;

/// @brief  Direct-form I 32.15
typedef struct
{
    int *b_k;
    int *a_k;
    int *y_k;
    int *x_k;
    size_t size_b;
    size_t size_a;
    size_t max_size;
} DirectForm1Q15; //

/// @brief  Direct-form II
typedef struct
{
    double *b_k;
    double *a_k;
    double *v_k;
    size_t size_b;
    size_t size_a;
    size_t max_size;
} DirectForm2;

/// @brief  Direct-form II 32.15
typedef struct
{
    int *b_k;
    int *a_k;
    int *v_k;
    size_t size_b;
    size_t size_a;
    size_t max_size;
} DirectForm2Q15;

/// @brief  Direct-form I Transposed
typedef struct
{
    double *b_k;
    double *a_k;
    double **s_a;
    double **s_b;
    size_t size_b;
    size_t size_a;
} DirectForm1Transposed;

/// @brief  Direct-form I Transposed 16.15
typedef struct
{
    int *b_k;
    int *a_k;
    int *s_a;
    int *s_b;
    size_t size_b;
    size_t size_a;
} DirectForm1TransposedQ15;

/// @brief  Direct-form II Transposed
typedef struct
{
    double *b_k;
    double *a_k;
    double *s_ab;
    size_t size_b;
    size_t size_a;
} DirectForm2Transposed;

/// @brief  Direct-form II Transposed 16.15
typedef struct
{
    int *b_k;
    int *a_k;
    int *s_ab;
    size_t size_b;
    size_t size_a;
} DirectForm2TransposedQ15;

void process_df1(DirectForm1 *ptr, double *x, double *y, size_t samples_num);
void process_df1_q15(DirectForm1Q15 *ptr, double *x, double *y, size_t samples_num);
void process_df2(DirectForm2 *ptr, double *x, double *y, size_t samples_num);
void process_df2_q15(DirectForm2Q15 *ptr, double *x, double *y, size_t samples_num);
void process_df1_transposed(DirectForm1Transposed *ptr, double *x, double *y, size_t samples_num);
void process_df1_transposed_q15(DirectForm1TransposedQ15 *ptr, double *x, double *y, size_t samples_num);
void process_df2_transposed(DirectForm2Transposed *ptr, double *x, double *y, size_t samples_num);
void process_df2_transposed_q15(DirectForm2TransposedQ15 *ptr, double *x, double *y, size_t samples_num);

DirectForm1 *create_df1(DigitalFilter *filter_ptr);
DirectForm1Q15 *create_df1_q15(DigitalFilter *filter_ptr);
DirectForm2 *create_df2(DigitalFilter *filter);
DirectForm2Q15 *create_df2_q15(DigitalFilter *filter_ptr);
DirectForm1Transposed *create_df1_transposed(DigitalFilter *filter);
DirectForm1TransposedQ15 *create_df1_transposed_q15(DigitalFilter *filter);
DirectForm2Transposed *create_df2_transposed(DigitalFilter *filter);

#endif // FILTERING_H
