#ifndef DIGITALFILTERS_H
#define DIGITALFILTERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AnalogFilters.h"
#include "Utils.h"

/**
 * @brief Structure to represent a digital filter.
 *
 * This structure is used to store the coefficients and parameters of
 * a digital filter, typically obtained after transforming an analog filter.
 * Coefficients in this filter are in form
 * Ref: Y [z](1 − 0.4z−1 + 0.2z−2) = X[z](0.2 + 0.4z−1 + 0.2z−2)
 * a0,a1,a2... = z-2,z-1,z0..., similarly for bk
 */
typedef struct
{
    double *a_k;           ///< Coefficients of the denominator polynomial.
    double *b_k;           ///< Coefficients of the numerator polynomial.
    int size_a;            ///< Number of coefficients in the denominator polynomial.
    int size_b;            ///< Number of coefficients in the numerator polynomial.
    int power_numerator;   ///< Power corresponding to the first vector coefficient
    int power_denominator; ///< Power corresponding to the first vector coefficient
} DigitalFilter;

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

/**
 * @brief Function calculating magnitude response of causal filter
 * M(\omega)=|H(e^{j\omega})| is a function of \omega that is called the magnitude response of the amplitude response of the system.
 * \begin{align*}M(\omega)=\sqrt{\big( \Re \{H(e^{j\omega}) \}\big)^{2}+\big( \Im \{H(e^{j\omega}) \}\big) ^{2}}\end{align*}
 * reference:
 * https://aleksandarhaber.com/magnitude-amplitude-and-phase-response-of-discrete-time-systems-and-filters/
 */
double *magnitude_response_digital_filter(DigitalFilter *p, double *magnitudes, int n, double fs);

/**
 * @brief Function calculating phase response of causal filter
 * phase response is then computed by solving this equation
 * (7) \begin{align*}\tan \left( \theta (\omega ) \right) = \frac{\Im \{H(e^{j\omega}) \}}{\Re \{H(e^{j\omega}) \}}\end{align*}
 * reference:
 * https://aleksandarhaber.com/magnitude-amplitude-and-phase-response-of-discrete-time-systems-and-filters/
 */
double *phase_response_digital_filter(DigitalFilter *p, double *phases, int n);

#endif // DIGITALFILTERS_H
