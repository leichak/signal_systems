#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BilinearTransform.h"
#include "DigitalFilters.h"

/**
 * @brief Frees the memory allocated for the digital filter structure.
 *
 * This function deallocates the memory used by the `DigitalFilter` structure,
 * including the memory for the coefficients of the numerator and denominator
 * polynomials.
 *
 * @param p A pointer to the `DigitalFilter` structure to be freed.
 */
void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

/**
 * @brief Calculates the magnitude response of a digital filter.
 *
 * To translate H(z) to analog response, we substitute ejω for z, where
 * ω is 2π*freq, with freq being the normalized frequency/sample rate:
 *
 * H(e^jω) = (a₀e^(0jω) + a₁e^(1jω) + ... + aₙe^(-kjω)) /
 *            (b₀e^(0jω) + b₁e^(1jω) + ... + bₖe^(-kjω))
 *
 * In this implementation, b represents the numerator and a the denominator.
 *
 * @param p A pointer to the `DigitalFilter` structure.
 * @param magnitudes An array to store the calculated magnitudes.
 * @param n The number of frequency points.
 * @param fs The sampling frequency.
 * @return A pointer to an array of frequency points.
 */
double *magnitude_response_digital_filter(DigitalFilter *p, double *magnitudes, int n, double fs)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }

    double T = 1 / fs; // Sampling period
    fill_n_with_step(w_k, n, -M_PI / fs, M_PI / fs);

    complex double numerator, denominator, hw = 0.0 + 0.0 * I;
    size_t size = max_int(p->size_a, p->size_b);

    for (size_t i = 0; i < n; i++) {
        numerator = 0.0 + 0.0 * I;
        denominator = 0.0 + 0.0 * I;

        for (size_t k = 0; k < size; k++) {
            if (k < p->size_b)
                numerator += p->b_k[k] * cexp(-w_k[i] * I * (double)(k + p->size_b - 1));
            if (k < p->size_a)
                denominator += p->a_k[k] * cexp(-w_k[i] * I * (double)(k + p->size_a - 1));
        }

        hw = numerator / denominator;
        magnitudes[i] = 20 * log10(cabs(hw)); // Convert to decibels
    }

    return w_k;
}

/**
 * @brief Calculates the phase response of a digital filter.
 *
 * The phase response is calculated using the equation:
 *
 * \[
 * \tan(\theta(\omega)) = \frac{\Im\{H(e^{j\omega})\}}{\Re\{H(e^{j\omega})\}}
 * \]
 *
 * @param p A pointer to the `DigitalFilter` structure.
 * @param phases An array to store the calculated phases.
 * @param n The number of frequency points.
 * @return A pointer to an array of frequency points.
 */
double *phase_response_digital_filter(DigitalFilter *p, double *phases, int n)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }

    fill_n_with_step(w_k, n, 0, 2 * M_PI); // Frequency points
    complex double numerator, denominator, hw = 0.0 + 0.0 * I;
    size_t size = max_int(p->size_a, p->size_b);

    for (size_t i = 0; i < n; i++) {
        numerator = 0.0 + 0.0 * I;
        denominator = 0.0 + 0.0 * I;

        for (size_t k = 0; k < size; k++) {
            if (k < p->size_b)
                numerator += p->b_k[k] * cexp(-(w_k[i]) * I * (double)(p->size_b - 1 - k));
            if (k < p->size_a)
                denominator += p->a_k[k] * cexp(-(w_k[i]) * I * (double)(p->size_a - 1 - k));
        }

        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        phases[i] = atan2(hw_im, hw_re); // Use atan2 for correct quadrant
    }

    return w_k;
}

/**
 * @brief Stability check for the filter.
 *
 * A filter is stable if its impulse response h_n decays to 0 as n approaches infinity.
 * A transfer function is stable if all its poles are inside the unit circle in the z-plane.
 *
 * Reference:
 * https://www.dsprelated.com/freebooks/filters/Stability_Revisited.html
 *
 * Stability can also be assessed using reflection coefficients, where:
 * All poles are inside the unit circle if all reflection coefficients are strictly between -1 and 1.
 *
 * Note: This section requires implementation based on your specific needs.
 */
