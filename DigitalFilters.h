#ifndef DIGITALFILTERS_H
#define DIGITALFILTERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Structure to represent a digital filter.
 *
 * This structure is used to store the coefficients and parameters of
 * a digital filter, typically obtained after transforming an analog filter.
 */
typedef struct
{
    double *a_k; ///< Coefficients of the denominator polynomial.
    double *b_k; ///< Coefficients of the numerator polynomial.
    double g_0;  ///< DC gain (gain at zero frequency).
    int size_a;  ///< Number of coefficients in the denominator polynomial.
    int size_b;  ///< Number of coefficients in the numerator polynomial.
} DigitalFilter;

/**
 * @brief Normalizes the digital filter coefficients based on the 0th numerator coefficient (b0).
 *
 * This function normalizes the coefficients of the digital filter by dividing
 * all coefficients by the value of b0 (the first coefficient of the numerator polynomial).
 *
 * @param p A pointer to the `DigitalFilter` structure to be normalized.
 *
 * Normalizing the coefficients in this way ensures that the filter's gain at DC
 * (or the low-frequency limit for low-pass filters) is properly adjusted.
 */
void normalize_to_b0(DigitalFilter *p);

/**
 * @brief Frees the memory allocated for the digital filter structure.
 *
 * This function deallocates the memory used by the `DigitalFilter` structure,
 * including the memory for the coefficients of the numerator and denominator
 * polynomials.
 *
 * @param p A pointer to the `DigitalFilter` structure to be freed.
 *
 * After calling this function, the pointer `p` should no longer be used
 * to refer to the filter without reinitializing it.
 */
void free_digital_filter(DigitalFilter *p);

#endif // DIGITALFILTERS_H
