#include <math.h>
#include <stdbool.h> // Added to use the 'bool' type
#include <stdio.h>
#include <stdlib.h>

#include "BilinearTransform.h"

// Function to perform the first step of the Horner method
// Reference: https://arxiv.org/pdf/2401.03071
void horner_step1_divide_sn_substitute(DigitalFilter *p, long double fs, long double prewarp_w0, bool prewarp)
{
    // Start with maximum powers of numerator and denominator
    int power_b_max = p->size_b - 1;
    int power_a_max = p->size_a - 1;

    // Determine the maximum power for transformation
    int power_max = (power_b_max > power_a_max) ? power_b_max : power_a_max;
    p->power_numerator -= power_max;   // Adjust numerator power
    p->power_denominator -= power_max; // Adjust denominator power

    // Take absolute values of powers
    p->power_numerator = abs(p->power_numerator);
    p->power_denominator = abs(p->power_denominator);

    // Allocate memory for new coefficients if needed
    if (p->power_denominator != p->size_a - 1) {
        long double *a_k = (long double *)calloc(p->power_denominator + 1, sizeof(long double));
        for (size_t k = 0; k < p->size_a; k++) {
            a_k[k] = p->a_k[k];
        }
        free(p->a_k);
        p->a_k = a_k;
        p->size_a = p->power_denominator + 1;
    }

    if (p->power_numerator != p->size_b - 1) {
        long double *b_k = (long double *)calloc(p->power_numerator + 1, sizeof(long double));
        for (size_t k = 0; k < p->size_b; k++) {
            b_k[k] = p->b_k[k];
        }
        free(p->b_k);
        p->b_k = b_k;
        p->size_b = p->power_numerator + 1;
    }

    long double K = prewarp ? prewarp_w0 / tan((prewarp_w0 / ((long double)2.0 * fs))) : (long double)2.0 * fs;

    // Scale numerator coefficients
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[k] /= pow(K, (long double)(p->size_b - 1 - k));
    }

    // Scale denominator coefficients
    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[k] /= pow(K, (long double)(p->size_a - 1 - k));
    }

    horner_step3_flip(p);
}

// Function to perform synthetic division to shift polynomial with n
// Reference: https://arxiv.org/pdf/2401.03071
void horner_step2_5_shift_polynomial_with_n(DigitalFilter *p, long double divisor)
{
    size_t size_b = p->size_b;
    size_t size_a = p->size_a;

    // Allocate memory for shifted coefficients
    long double *division_b = (long double *)calloc(size_b, sizeof(long double));
    long double *division_a = (long double *)calloc(size_a, sizeof(long double));

    // Flip coefficients for processing
    for (size_t k = 0; k < size_b; k++) {
        division_b[k] = p->b_k[size_b - 1 - k];
    }
    for (size_t k = 0; k < size_a; k++) {
        division_a[k] = p->a_k[size_a - 1 - k];
    }

    // Shift coefficients for numerator
    size_t up_to = size_b;
    while (up_to > 0) {
        for (size_t i = 1; i < up_to; i++) {
            division_b[i] += division_b[i - 1] * divisor;
        }
        up_to--;
    }

    // Shift coefficients for denominator
    up_to = size_a;
    while (up_to > 0) {
        for (size_t i = 1; i < up_to; i++) {
            division_a[i] += division_a[i - 1] * divisor;
        }
        up_to--;
    }

    // Write back in original order
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[size_b - 1 - k] = division_b[k];
    }
    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[size_a - 1 - k] = division_a[k];
    }

    free(division_b);
    free(division_a);
}

// Function to flip coefficients
void horner_step3_flip(DigitalFilter *p)
{
    // Flip numerator coefficients
    long double *p1 = p->b_k;
    long double *p2 = p->b_k + p->size_b - 1;

    while (p1 < p2) {
        long double t = *p1;
        *p1 = *p2;
        *p2 = t;

        p1++;
        p2--;
    }

    // Flip denominator coefficients
    p1 = p->a_k;
    p2 = p->a_k + p->size_a - 1;

    while (p1 < p2) {
        long double t = *p1;
        *p1 = *p2;
        *p2 = t;

        p1++;
        p2--;
    }
}

// Function to scale polynomial zeros by 2
// Reference: https://arxiv.org/pdf/2401.03071
void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p)
{
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[k] *= pow(0.5, (long double)k);
    }

    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[k] *= pow(0.5, (long double)k);
    }
}

// Function to normalize filter coefficients to make the system causal
void horner_step6_make_causal_normalize_to_a0(DigitalFilter *p, long double fs)
{
    long double norm = p->a_k[p->size_a - 1];

    // Normalize coefficients
    for (size_t k = 0; k < p->size_b; k++)
        p->b_k[k] *= ((long double)1.0 / norm);
    for (size_t k = 0; k < p->size_a; k++)
        p->a_k[k] *= ((long double)1.0 / norm);
}

// Main function to perform bilinear transformation using the Horner method
DigitalFilter *bilinear_transform_horner_method(AnalogFilter *p, double fs, double w0)
{
    DigitalFilter *p_d = (DigitalFilter *)malloc(sizeof(DigitalFilter));
    if (p_d == NULL)
        return NULL;

    // Allocate memory for coefficients
    p_d->a_k = (long double *)calloc(p->size_a, sizeof(long double));
    p_d->b_k = (long double *)calloc(p->size_b, sizeof(long double));

    // Copy coefficients from analog filter
    for (size_t k = 0; k < p->size_b; k++)
        p_d->b_k[k] = p->b_k[k];
    for (size_t k = 0; k < p->size_a; k++)
        p_d->a_k[k] = p->a_k[k];

    // Initialize sizes and powers
    p_d->size_a = p->size_a;
    p_d->size_b = p->size_b;
    p_d->power_numerator = 0;
    p_d->power_denominator = 0;

    // Apply Horner method steps
    horner_step1_divide_sn_substitute(p_d, (long double)fs, (long double)w0, true);
    horner_step2_5_shift_polynomial_with_n(p_d, 1.0);
    horner_step3_flip(p_d);
    horner_step4_scale_polynomial_zeros_by_2(p_d);
    horner_step2_5_shift_polynomial_with_n(p_d, -1.0);
    horner_step6_make_causal_normalize_to_a0(p_d, 1.0);

    return p_d;
}

// Function to convert analog frequency to digital frequency
double wa_2_wd(double wa, double T)
{
    return (2.0 / T) * atan(wa * (T / 2.0));
}

// Function to check the stability of the system
int stabilitycheck(double *A, int length)
{
    int N = length - 1; // Order of A(z)
    int stable = 1;     // Assume stable unless shown otherwise

    // Check for stability by examining reflection coefficients
    for (int i = N - 1; i > -1; i--) {
        double rci = A[i + 1]; // Get reflection coefficient

        if (fabs(rci) >= 1.0) { // If coefficient is outside [-1, 1], unstable
            stable = 0;
            break;
        }
    }
    return stable; // Return stability status
}

// Function to convert double frequency to radians
double double_frequency_2_radians(double f)
{
    return 2.0 * M_PI * f; // Convert frequency to radians
}

// Function to check frequency response
void check_frequency_response(DigitalFilter *p, double f)
{
    // Frequency response computation (to be implemented)
    // Placeholder for frequency response check
}

// Add any additional functions needed for your filter processing
