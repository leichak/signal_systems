#include "AnalogFilters.h"
#include "DigitalFilters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *transform_polynomial(int Order, double beta, double gamma, double delta, double alpha, double *coefficients_original)
{
    double *beta_powers = malloc((Order + 1) * sizeof(double));
    double *gamma_powers = malloc((Order + 1) * sizeof(double));
    double *transformed_coefficients = malloc((Order + 1) * sizeof(double));
    double *coefficients = malloc((Order + 1) * sizeof(double));
    memcpy(coefficients, coefficients_original, sizeof(double) * (Order + 1));

    if (!beta_powers || !gamma_powers || !transformed_coefficients) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(beta_powers);
        free(gamma_powers);
        free(transformed_coefficients);
        return NULL;
    }

    beta_powers[0] = gamma_powers[0] = 1.0;

    for (int j = 1; j <= Order; j++) {
        beta_powers[j] = beta * beta_powers[j - 1];
        gamma_powers[j] = gamma * gamma_powers[j - 1];
    }

    for (int k = 0; k < Order; k++) {
        double sum = 0.0;
        for (int j = 0; j <= Order - k; j++) {
            coefficients[j] = ((Order - k - j) * delta * coefficients[j] + (j + 1) * alpha * coefficients[j + 1]) / (k + 1);
            sum += beta_powers[j] * gamma_powers[Order - k - j] * coefficients[j];
        }
        transformed_coefficients[Order - k] = sum;
    }

    transformed_coefficients[0] = coefficients[0];

    free(beta_powers);
    free(gamma_powers);

    return transformed_coefficients;
}

DigitalFilter *transform_analog_to_digital(AnalogFilter *p_a)
{
    DigitalFilter *p_d = malloc(sizeof(DigitalFilter));
    if (!p_d) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    p_d->b_k = transform_polynomial(p_a->order_numerator, -1.0, 1.0, 1.0, 1.0, p_a->b_k);
    p_d->a_k = transform_polynomial(p_a->order_denominator, -1.0, 1.0, 1.0, 1.0, p_a->a_k);

    p_d->size_a = p_a->size_a;
    p_d->size_b = p_a->size_b;

    if (!p_d->b_k || !p_d->a_k) {
        free(p_d);
        fprintf(stderr, "Polynomial transformation failed.\n");
        return NULL;
    }

    return p_d;
}

// // https://arxiv.org/pdf/2401.03071
void horner_step1_divide_sn_substitute(DigitalFilter *p)
{
    // divide by the highest order and substitute s for 2fl/x
    // input is  b x^i ... i=start+0..start+1...start+2...len()
    // size will stay the same so start power is sufficient?
    int power_b = (p->power_numerator + (p->size_b - 1));
    int power_a = (p->power_denominator + (p->size_a - 1));

    // Divide by a power
    int max_power = power_b > power_a ? power_b : power_a;
    p->power_numerator -= max_power;
    p->power_denominator -= max_power;
}

// https://arxiv.org/pdf/2401.03071
void horner_decrease_by_n_with_synthetic_division(DigitalFilter *p, int n)
{

    // Next we can decrease the zeros by 1 using synthetic division
}

void horner_step3_replace_all_zeros_with_reciprocals(DigitalFilter *p)
{
    // just flipping coefficients
    double *p1 = p->b_k;
    double *p2 = p->b_k + p->size_b;

    while (p1 != p2) {
        double t = *p1;
        *p1 = *p2;
        *p2 = t;

        p1++;
        p2--;
    }

    // just flipping coefficients
    p1 = p->a_k;
    p2 = p->a_k + p->size_a;

    while (p1 != p2) {
        double t = *p1;
        *p1 = *p2;
        *p2 = t;

        p1++;
        p2--;
    }
}

// // https://arxiv.org/pdf/2401.03071
// void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p)
// {
//     /*
//     Scale the polynomial zeroes by 2. Note can you could scale by
//         either 1
//     2 or 2 for the orders of power since the polynomials are in both the numerator and
//     denominator.
//     */
// }

// // https://arxiv.org/pdf/2401.03071
// void horner_step5_increase_by_1_with_synthetic_division(DigitalFilter *p)
// {
//     //  Increase all polynomial zeros by 1 using synthetic division
// }

// // https://arxiv.org/pdf/2401.03071
// void billinear_transform_horner_method()
// {
// }
