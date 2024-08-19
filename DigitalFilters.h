#ifndef DIGITALFILTERS_H
#define DIGITALFILTERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Structure to represent digital filter.
 */
typedef struct
{
    double *a_k; ///< Coefficients of the denominator polynomial.
    double *b_k; ///< Coefficients of the numerator polynomial.
    double g_0;  ///< DC gain.
    int size_a;  ///< Number of coefficients in the denominator.
    int size_b;  ///< Number of coefficients in the numerator.
} DigitalFilter;

/** @brief This function norms digital filter to b0, a.k.a multiplies
 * coefficients by 1/b0
 *
 */
void normalize_to_b0(DigitalFilter *p);

void free_digital_filter(DigitalFilter *p);

#endif
