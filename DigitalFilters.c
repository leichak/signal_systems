#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BilinearTransform.h"
#include "DigitalFilters.h"

void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

/* reference: https://www.earlevel.com/main/2016/12/01/evaluating-filter-frequency-response/
 *
 *  To translate H(z) to analog response, we substitute ejw for z, where w is 2pi*freq,
 *  with freq being the normalized frequency/sample rate thus:
 *
 *  H(ejw) = (a0e^0jw + a1e^1jw + ... + ane^-kjw) / (abe^0jw + abe^1jw + ... + abe^-kjw) /
 *  (in our implementation b is numerator and a denominator)
 *
 */

double *
magnitude_response_digital_filter(DigitalFilter *p, double *magnitudes, int n, double fs)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }
    double T = 1 / fs;
    fill_n_with_step(w_k, n, -M_PI / fs, M_PI / fs);

    complex double numerator, denominator, hw = 0.0 + 0.0 * I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        numerator = 0 + 0 * I;
        denominator = 0 + 0 * I;
        for (size_t k = 0; k < size; k++) {
            if (k < p->size_b)
                numerator += p->b_k[k] * cexp(-w_k[i] * I * (double)(k + p->size_b - 1));
            if (k < p->size_a)
                denominator += p->a_k[k] * cexp(-w_k[i] * I * (double)(k + p->size_a - 1));
        }
        hw = numerator / denominator;
        magnitudes[i] = 20 * log10(cabs(hw));
    }

    return w_k;
}

double *phase_response_digital_filter(DigitalFilter *p, double *phases, int n)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }
    fill_n_with_step(w_k, n, 0, 2 * M_PI);
    for (size_t k = 0; k < n; k++) {
        // w_k[k] = rad_s_2_hz(wa_2_wd(w_k[k], 1 / fs)); // frequency warping compensation + conversion to Hz
        // w_k[k] = rad_s_2_hz(w_k[k]);
    }
    complex double numerator, denominator, hw = 0 * 0.0I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        numerator = 0.0;
        denominator = 0.0;
        for (size_t k = 0; k < size; k++) {
            if (k < p->size_b)
                numerator += p->b_k[k] * cexp(-(w_k[i]) * I * (double)(p->size_b - 1 - k));
            if (k < p->size_a)
                denominator -= p->a_k[k] * cexp(-(w_k[i]) * I * (double)(p->size_a - 1 - k));
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        phases[i] = atan(hw_im / hw_re);
    }

    return w_k;
}

/** https://www.dsprelated.com/freebooks/filters/Stability_Revisited.html
 * As defined earlier filter is said to be stable if its impulse response hn decays to 0 as
 * n goes to infinity. In the terms of poles and zeros, an irreducible filter transfer function
 * is stable if and only if all its poles are inside the unit circle in the z plan. This is
 * because transfer function is the z transform of the impulse response, and if there is observable
 * pole outside the unit circle then there is an exponentially increasing component of the
 * impulse response. So all poles need to be inside the unit circle.
 *
 * Poles on the unit circle may be called marginally stable. The impulse response component
 * corresponding to a single pole on the unit circle never decays but neither does it grow.
 *
 * Apart from checking that all poles are residing inside unit circle there is faster method
 * based on filter reflection coefficients. It is mathematical fact that all poles of a recursive
 * filter are inside the unit curcle if and only if all its reflection coefficients are strictly
 * between -1 and 1.
 *
 * Step-down procedure, Schur Cohn stability test or Durbin recursion (Levinson algorithm)
 */
