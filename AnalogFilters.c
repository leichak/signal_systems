#include "AnalogFilters.h"

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
void butterworth_coefficients(int order, double *coefficients)
{
    double gamma = M_PI / (2.0 * order);

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
    f->b_k = (double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL) {
        free(f);
        return NULL;
    }
    f->b_k[0] = 1.0;

    f->size_a = n + 1;
    f->a_k = (double *)malloc(f->size_a * sizeof(double));
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
        p->a_k[i] *= pow(1.0 / wc, i);
    }
}

// Performs a high-pass frequency transformation on the filter
void transform_to_high_pass(AnalogFilter *p, double wc)
{
    int order_numerator = p->order_numerator;
    int order_denominator = p->order_denominator;

    // Scale the coefficients by wc
    for (size_t i = 1; i < p->size_a; i++) {
        p->a_k[i] *= pow(wc, i);
    }
    for (size_t i = 1; i < p->size_b; i++) {
        p->b_k[i] *= pow(wc, i);
    }

    // Adjust order to account for the high-pass transformation
    p->order_denominator += order_numerator;
    p->size_a = p->order_denominator + 1;
    p->order_numerator += order_denominator;
    p->size_b = p->order_numerator + 1;

    double *b = (double *)malloc(p->size_b * sizeof(double));
    double *a = (double *)malloc(p->size_a * sizeof(double));
    size_t size = max_int(p->size_a, p->size_b);

    // Perform order adjustments (shifting polynomial coefficients)
    for (size_t i = 0; i < size; i++) {
        if (i < p->size_b) {
            b[i] = (i >= order_denominator) ? p->b_k[i - order_denominator] : 0.0;
        }
        if (i < p->size_a) {
            a[i] = (i >= order_numerator) ? p->a_k[i - order_numerator] : 0.0;
        }
    }

    // Reverse the orders (required for high-pass transformation)
    for (int i = 0, j = p->size_b - 1; i < j && i != j; i++, j--) {
        double t = b[j];
        b[j] = b[i];
        b[i] = t;
    }
    for (int i = 0, j = p->size_a - 1; i < j && i != j; i++, j--) {
        double t = a[j];
        a[j] = a[i];
        a[i] = t;
    }

    free(p->a_k);
    free(p->b_k);
    p->a_k = a;
    p->b_k = b;
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
                denominator += p->a_k[k] * cpow((w_k[i] / M_PI) * I, (double)(k));
            if (k < p->size_b)
                numerator += p->b_k[k] * cpow((w_k[i] / M_PI) * I, (double)(k));
        }
        hw = numerator / denominator;
        double hw_re = creal(hw);
        double hw_im = cimag(hw);
        magnitudes[i] = sqrt(pow(hw_re, 2) + pow(hw_im, 2));
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
