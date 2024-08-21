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

void free_causal_digital_filter(DigitalFilterCausal *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

void normalize_to_b0(DigitalFilter *p)
{
    if (p == NULL) {
        return;
    }

    double a_0 = p->a_k[0];

    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] /= a_0;
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] /= a_0;
    }
}

void normalize_to_b0_causal(DigitalFilterCausal *p)
{
    if (p == NULL) {
        return;
    }

    double a_0 = p->a_k[p->size_a - 1];

    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] /= a_0;
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] /= a_0;
    }
}

DigitalFilterCausal *make_causal(DigitalFilter *p)
{
    if (!p) {
        return NULL;
    }

    DigitalFilterCausal *p_d_c = (DigitalFilterCausal *)malloc(sizeof(DigitalFilterCausal));
    assert(p_d_c != NULL);

    size_t max_size = p->size_a > p->size_b ? p->size_a : p->size_b;
    p_d_c->a_k = (double *)calloc(max_size, sizeof(double));
    p_d_c->b_k = (double *)calloc(max_size, sizeof(double));
    assert(p_d_c->a_k != NULL && p_d_c->b_k != NULL);

    memcpy(p_d_c->a_k, p->a_k, sizeof(double) * p->size_a);
    memcpy(p_d_c->b_k, p->b_k, sizeof(double) * p->size_b);

    p_d_c->size_a = max_size;
    p_d_c->size_b = max_size;

    return p_d_c;
}

int frequency_response_causal_digital_filter(DigitalFilterCausal *p, double *magnitudes, int n)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return -1;
    }
    fill_n_with_step(w_k, n, -M_PI, M_PI);
    complex numerator, denominator, hw = 0 * 0.0I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < size; i++) { // just calculate
            if (k < p->size_a)
                denominator += p->a_k[p->size_a - 1 - k] * cexp(-k * w_k[i]);
            if (k < p->size_b)
                numerator += p->b_k[p->size_b - 1 - k] * cexp(-k * w_k[i]);
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        magnitudes[i] = sqrt(powf(hw_re, 2) + powf(hw_im, 2));
    }

    free(w_k);

    return 0;
}

int phase_response_causal_digital_filter(DigitalFilterCausal *p, double *phases, int n)
{
    double *w_k = (double *)calloc(n, sizeof(double));
    if (w_k == NULL) {
        return -1;
    }
    fill_n_with_step(w_k, n, -M_PI, M_PI);
    complex numerator, denominator, hw = 0 * 0.0I;
    size_t size = max_int(p->size_a, p->size_b);
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < size; i++) { // just calculate
            if (k < p->size_a)
                denominator += p->a_k[p->size_a - 1 - k] * cexp(-k * w_k[i]);
            if (k < p->size_b)
                numerator += p->b_k[p->size_b - 1 - k] * cexp(-k * w_k[i]);
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        phases[i] = atan(hw_im / hw_re);
    }

    free(w_k);

    return 0;
}
