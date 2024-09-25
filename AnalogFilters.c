#include "AnalogFilters.h"
#include "Utils.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Computes the Butterworth filter gain for a low-pass filter
float butter_gain_low_pass(int n, float g0, float wc, float w)
{
    return pow(g0, 2) / (1 + pow((w / wc), 2 * n));
}

// Computes the logarithmic Butterworth filter gain for a low-pass filter
float butter_gain_low_pass_log10(int n, float g0, float wc, float w)
{
    return 10.0f * log10(pow(g0, 2) / (1 + pow((w / wc), 2 * n)));
}

// Frees allocated memory for the coefficients of the analog filter
void free_analog_filter(AnalogFilter *p)
{
    if (p != NULL) {
        // Free the coefficients of the denominator polynomial (a_k)
        if (p->a_k != NULL) {
            free(p->a_k);
            p->a_k = NULL;
        }

        // Free the coefficients of the numerator polynomial (b_k)
        if (p->b_k != NULL) {
            free(p->b_k);
            p->b_k = NULL;
        }

        // Free the AnalogFilter structure itself
        free(p);
    }
}

// Generates the coefficients for a Butterworth filter of a given order
void butterworth_coefficients(int order, long double *coefficients)
{
    long double gamma = M_PI / ((long double)2.0 * (long double)order);

    // Initialize the first coefficient (a0)
    coefficients[0] = 1.0;

    // Compute remaining coefficients using the product formula
    for (int k = 1; k <= order; ++k) {
        coefficients[k] = coefficients[k - 1] * (cos((k - 1) * gamma) / sin(k * gamma));
    }

    // Exploit symmetry: a_k = a_(n-k)
    for (int k = 0; k <= order / 2; ++k) {
        coefficients[order - k] = coefficients[k];
    }
}

// Generates an analog filter based on the given specifications
AnalogFilter *generate_analog_filter(int n, double wc, FilterTypes filter_type, BandType band_type)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    switch (filter_type) {
    case BUTTERWORTH:
        assert(generate_butterworth(n, f) != NULL);
        break;
    }

    // Apply the appropriate frequency transformation
    switch (band_type) {
    case LOWPASS:
        transform_to_low_pass(f, wc);
        break;
    case HIGHPASS:
        transform_to_high_pass(f, wc);
        break;
    default:
        transform_to_low_pass(f, wc);
    }

    return f;
}

// Generates a Butterworth filter and initializes its coefficients
void *generate_butterworth(int n, AnalogFilter *f)
{
    f->g_0 = 1.0;
    f->order_denominator = n;
    f->order_numerator = 0;
    f->size_b = f->order_numerator + 1;
    f->b_k = (long double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL) {
        free(f);
        return NULL;
    }
    f->b_k[0] = 1.0;

    f->size_a = n + 1;
    f->a_k = (long double *)malloc(f->size_a * sizeof(double));
    if (f->a_k == NULL) {
        free(f->b_k);
        free(f);
        return NULL;
    }
    butterworth_coefficients(f->order_denominator, f->a_k);

    return f; // Return a non-null pointer to indicate successful allocation
}

// Performs a low-pass frequency transformation on the filter
void transform_to_low_pass(AnalogFilter *p, double wc)
{
    for (int i = 1; i < p->size_a; i++) {
        p->a_k[i] *= pow((long double)1.0 / (long double)wc, (long double)i);
    }
}

// Performs a high-pass frequency transformation on the filter
void transform_to_high_pass(AnalogFilter *p, double wc)
{
}

// Normalizes the filter coefficients to their maximum
void normalize_to_max(AnalogFilter *p)
{
    double a_max = p->a_k[0];
    double b_max = p->b_k[0];

    for (int i = 1; i < p->size_b; i++) {
        if (p->b_k[i] > b_max)
            b_max = p->b_k[i];
    }

    for (int i = 1; i < p->size_a; i++) {
        if (p->a_k[i] > a_max)
            a_max = p->a_k[i];
    }

    p->g_0 = a_max / b_max;

    for (int i = 1; i < p->size_a; i++) {
        p->a_k[i] /= p->g_0;
    }
}

double *magnitude_response_analog_filter(AnalogFilter *p, double *magnitudes, int n)
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
        for (size_t k = 0; k < size; k++) {
            if (k < p->size_a)
                denominator += p->a_k[k] * cpow((w_k[i]) * I, (double)(k));
            if (k < p->size_b)
                numerator += p->b_k[k] * cpow((w_k[i]) * I, (double)(k));
        }
        hw = numerator / denominator;
        magnitudes[i] = 20 * log10(cabs(hw));
    }

    return w_k;
}

double *phase_response_analog_filter(AnalogFilter *p, double *phases, int n)
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
        for (size_t k = 0; k < size; k++) {
            if (k < p->size_a)
                denominator += p->a_k[k] * cpow((w_k[i] / M_PI) * I, (double)(k));
            if (k < p->size_b)
                numerator += p->b_k[k] * cpow((w_k[i] / M_PI) * I, (double)(k));
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        phases[i] = atan(hw_im / hw_re);
    }

    return w_k;
}
