#ifndef DIGITALFILTERS_H
#define DIGITALFILTERS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AnalogFilters.h"
#include "Utils.h"

/**
 * @brief Structure representing a digital filter.
 *
 * This structure stores the coefficients and parameters of a digital filter,
 * typically obtained after transforming an analog filter. The coefficients are
 * represented in the form:
 *
 * Y[z](1 − 0.4z⁻¹ + 0.2z⁻²) = X[z](0.2 + 0.4z⁻¹ + 0.2z⁻²)
 *
 * where a₀, a₁, a₂... correspond to z⁻², z⁻¹, z⁰..., and similarly for bₖ.
 */
typedef struct
{
    long double *a_k;      ///< Coefficients of the denominator polynomial.
    long double *b_k;      ///< Coefficients of the numerator polynomial.
    int size_a;            ///< Number of coefficients in the denominator polynomial.
    int size_b;            ///< Number of coefficients in the numerator polynomial.
    int power_numerator;   ///< Maximum power in the numerator polynomial.
    int power_denominator; ///< Maximum power in the denominator polynomial.
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
 * After calling this function, the pointer `p` should not be used
 * to refer to the filter without reinitializing it.
 */
void free_digital_filter(DigitalFilter *p);

/**
 * @brief Calculates the magnitude response of a causal filter.
 *
 * The magnitude response is defined as:
 *
 * \[
 * M(\omega) = |H(e^{j\omega})| = \sqrt{(\Re\{H(e^{j\omega})\})^2 + (\Im\{H(e^{j\omega})\})^2}
 * \]
 *
 * @param p A pointer to the `DigitalFilter` structure.
 * @param magnitudes An array to store the calculated magnitudes.
 * @param n The number of frequency points.
 * @param fs The sampling frequency.
 * @return A pointer to the array containing the magnitude response.
 *
 * Reference:
 * https://aleksandarhaber.com/magnitude-amplitude-and-phase-response-of-discrete-time-systems-and-filters/
 */
double *magnitude_response_digital_filter(DigitalFilter *p, double *magnitudes, int n, double fs);

/**
 * @brief Calculates the phase response of a causal filter.
 *
 * The phase response is computed using the equation:
 *
 * \[
 * \tan(\theta(\omega)) = \frac{\Im\{H(e^{j\omega})\}}{\Re\{H(e^{j\omega})\}}
 * \]
 *
 * @param p A pointer to the `DigitalFilter` structure.
 * @param phases An array to store the calculated phases.
 * @param n The number of frequency points.
 * @return A pointer to the array containing the phase response.
 *
 * Reference:
 * https://aleksandarhaber.com/magnitude-amplitude-and-phase-response-of-discrete-time-systems-and-filters/
 */
double *phase_response_digital_filter(DigitalFilter *p, double *phases, int n);

#endif // DIGITALFILTERS_H
