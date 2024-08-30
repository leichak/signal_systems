#include "DigitalFilters.h"
#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "Utils.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

double *magnitude_response_digital_filter(DigitalFilter *p, double *magnitudes, int n, double fs)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }
    fill_n_with_step(w_k, n, -M_PI, M_PI);
    complex double numerator, denominator, hw = 0.0 + 0.0 * I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        numerator = 0 + 0 * I;
        denominator = 0 + 0 * I;
        for (size_t k = 0; k < size; k++) {
            if (k < p->size_b)
                numerator += p->b_k[k] * cexp(-(w_k[i]) * I * (double)(p->power_numerator - k));
            if (k < p->size_a)
                denominator += p->a_k[k] * cexp(-(10 * w_k[i]) * I * (double)(p->power_denominator - k));
        }
        hw = numerator / denominator;
        magnitudes[i] = cabs(hw);
    }

    return w_k;
}

double *phase_response_digital_filter(DigitalFilter *p, double *phases, int n)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return NULL;
    }
    fill_n_with_step(w_k, n, -M_PI, M_PI);
    complex double numerator, denominator, hw = 0 * 0.0I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        numerator = 0.0;
        denominator = 0.0;
        for (size_t k = 0; k < size; k++) { // just calculate
            if (k < p->size_a)
                denominator += p->a_k[p->size_a - 1 - k] * cexp(-k * w_k[i] * I);
            if (k < p->size_b)
                numerator += p->b_k[p->size_b - 1 - k] * cexp(-k * w_k[i] * I);
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        phases[i] = atan(hw_im / hw_re);
    }

    return w_k;
}
