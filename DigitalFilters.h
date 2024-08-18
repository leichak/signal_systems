#ifndef DIGITALFILTERS_H
#define DIGITALFILTERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Structure to represent digital filter.
 */
typedef struct {
  double *a_k;  ///< Coefficients of the denominator polynomial.
  double *b_k;  ///< Coefficients of the numerator polynomial.
  double g_0;   ///< DC gain.
  int size_a;   ///< Number of coefficients in the denominator.
  int size_b;   ///< Number of coefficients in the numerator.
} DigitalFilter;

void free_digital_filter(DigitalFilter *p);

#endif
