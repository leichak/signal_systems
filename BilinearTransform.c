#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "BilinearTransform.h"

// // https://arxiv.org/pdf/2401.03071
void horner_step1_divide_sn_substitute(DigitalFilter *p, double fs, double prewarp_w0, bool prewarp)
{
    // coefficients here should be in increasing order, 0,1,2,3,4 power
    // dividing by the 1/(s^ max_power) makes 4rth power 1 and, 0 power -4
    // substituting 2fl / x then will essentially flip coefficients and leave
    // multiplied by (1/2fl)^power, so powers are 4,3,2,1,0,
    // we should multiply and flip again to make it ok

    // start max power
    int power_b_max = p->size_b - 1; // start power
    int power_a_max = p->size_a - 1; // start power,

    // powers vars will serve temp for some transformations, later they should reflect min power
    // multiplying by 1/s
    int power_max = power_b_max > power_a_max ? power_b_max : power_a_max;
    p->power_numerator -= power_max;   // s0 -> s-2
    p->power_denominator -= power_max; // s0 -> s-2, s2 => s0

    // substituting (2fs/x)
    p->power_numerator = abs(p->power_numerator);
    p->power_denominator = abs(p->power_denominator);

    // we can assume that if there is difference with final power and size we need to allocate
    if (p->power_denominator != p->size_a - 1) {
        double *a_k = (double *)calloc(p->power_denominator + 1, sizeof(double));
        for (size_t k = 0; k < p->size_a; k++) {
            a_k[k] = p->a_k[k];
        }
        free(p->a_k);
        p->a_k = a_k;
        p->size_a = p->power_denominator + 1;
    }

    if (p->power_numerator != p->size_b - 1) {
        double *b_k = (double *)calloc(p->power_numerator + 1, sizeof(double));
        for (size_t k = 0; k < p->size_b; k++) {
            b_k[k] = p->b_k[k];
        }
        free(p->b_k);
        p->b_k = b_k;
        p->size_b = p->power_numerator + 1;
    }

    double K;

    if (prewarp)
        K = prewarp_w0 / tan((prewarp_w0 / (2.0 * fs)));
    else
        K = 2.0 * fs;

    // fl is sampling frequency, multiplication
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[k] /= powf(K, ((long double)p->size_b - 1 - (long double)k));
    }

    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[k] /= powf(K, ((long double)p->size_a - 1 - (long double)k));
    }

    horner_step3_flip(p);
}

// https://arxiv.org/pdf/2401.03071
void horner_step2_5_shift_polynomial_with_n(DigitalFilter *p, long double divisor)
{
    // Next we can decrease the zeros by 1 using synthetic division
    // We need to take remainders
    // Assume that power are increasing with indexes (0,1,2,3,4), that is intuitive
    size_t size_b = p->size_b;
    size_t size_a = p->size_a;
    double *division_b = (double *)calloc(size_b, sizeof(double));
    double *division_a = (double *)calloc(size_a, sizeof(double));
    for (size_t k = 0; k < size_b; k++) { // flip them for better code
        division_b[k] = p->b_k[size_b - 1 - k];
    }
    for (size_t k = 0; k < size_a; k++) { // flip them for better code
        division_a[k] = p->a_k[size_a - 1 - k];
    }

    size_t up_to = size_b;
    while (up_to > 0) {
        for (size_t i = 1; i < up_to; i++) {
            // division
            division_b[i] += division_b[i - 1] * divisor;
        }
        up_to -= 1;
    }

    up_to = size_a;
    while (up_to > 0) {

        for (size_t i = 1; i < up_to; i++) {
            // division
            division_a[i] += division_a[i - 1] * divisor;
        }
        up_to -= 1;
    }

    // write in original order // flip again // powers 0,1,2,3,4,5
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[size_b - 1 - k] = division_b[k];
    }
    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[size_a - 1 - k] = division_a[k];
    }

    free(division_b);
    free(division_a);
}

void horner_step3_flip(DigitalFilter *p)
{
    // just flipping coefficients
    long double *p1 = p->b_k;
    long double *p2 = p->b_k + p->size_b - 1;

    while (p1 < p2) {
        long double t = *p1;
        *p1 = *p2;
        *p2 = t;

        p1++;
        p2--;
    }

    // just flipping coefficients
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

// https://arxiv.org/pdf/2401.03071
void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p)
{
    for (size_t k = 0; k < p->size_b; k++) {
        p->b_k[k] *= powf(0.5, (long double)k);
    }

    for (size_t k = 0; k < p->size_a; k++) {
        p->a_k[k] *= powf(0.5, (long double)k);
    }
}

void horner_step6_make_causal_normalize_to_a0(DigitalFilter *p, double fs)
{
    // not sure whether i need it
    int power_b_max = p->size_b - 1;
    int power_a_max = p->size_a - 1;

    // multiplying by 1/s^n max
    int power_max = power_b_max > power_a_max ? power_b_max : power_a_max;
    // p->power_numerator -= power_max;
    // p->power_denominator -= power_max;

    p->power_numerator = p->power_numerator - 2 * (power_max);
    p->power_denominator = p->power_denominator - 2 * (power_max);

    // for (size_t k = 0; k < p->size_a; k++) {
    //     p->a_k[k] *= 1.0 / powf(fs, -k);
    // }

    // for (size_t k = 0; k < p->size_b; k++) {
    //     p->b_k[k] *= 1.0 / powf(fs, -k);
    // }

    // reflect sample rate

    // make y0 to be 1
    double norm = p->a_k[p->size_a - 1];
    for (size_t k = 0; k < p->size_b; k++)
        p->b_k[k] *= (1.0 / norm);
    for (size_t k = 0; k < p->size_a; k++)
        p->a_k[k] *= (1.0 / norm);
}

DigitalFilter *bilinear_transform_horner_method(AnalogFilter *p, double fs, double w0)
{
    DigitalFilter *p_d = (DigitalFilter *)malloc(sizeof(DigitalFilter));
    if (p_d == NULL)
        return NULL;
    p_d->a_k = (double *)calloc(p->size_a, sizeof(double));
    p_d->b_k = (double *)calloc(p->size_b, sizeof(double));
    for (size_t k = 0; k < p->size_b; k++)
        p_d->b_k[k] = p->b_k[k];
    for (size_t k = 0; k < p->size_a; k++)
        p_d->a_k[k] = p->a_k[k];
    p_d->size_a = p->size_a;
    p_d->size_b = p->size_b;
    p_d->power_numerator = 0;
    p_d->power_denominator = 0;

    // Input coefficients from analog are in order 0,1...n power
    horner_step1_divide_sn_substitute(p_d, fs, w0, true);
    horner_step2_5_shift_polynomial_with_n(p_d, 1.0);
    horner_step3_flip(p_d);
    horner_step4_scale_polynomial_zeros_by_2(p_d);
    horner_step2_5_shift_polynomial_with_n(p_d, -1);
    horner_step6_make_causal_normalize_to_a0(p_d, 1.0);

    return p_d;
}

double wa_2_wd(double wa, double T)
{
    return (2.0 / T) * atan(wa * (T / 2.0));
}

int stabilitycheck(double *A, int length)
{
    int N = length - 1; // Order of A(z)
    int stable = 1;     // Assume stable unless shown otherwise

    // Ensure A is a column vector (already handled as a pointer in C)
    for (int i = N - 1; i >= 0; i--) {
        double rci = A[i + 1]; // Get reflection coefficient

        if (fabs(rci) >= 1) {
            stable = 0;
            return stable; // System is unstable
        }

        // Update A array with new coefficients
        for (int j = 0; j < i; j++) {
            A[j] = (A[j] - rci * A[i + 1 - j]) / (1 - rci * rci);
        }
    }

    return stable; // Return 1 if stable, 0 if unstable
}
