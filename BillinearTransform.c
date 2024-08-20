#include "AnalogFilters.h"
#include "DigitalFilters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *transform_polynomial(int Order, double beta, double gamma, double delta, double alpha, double *coefficients)
{
    double *beta_powers = malloc((Order + 1) * sizeof(double));
    double *gamma_powers = malloc((Order + 1) * sizeof(double));
    double *transformed_coefficients = malloc((Order + 1) * sizeof(double));

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

DigitalFilter *transform_analog_to_digital(AnalogFilter *pa)
{
    DigitalFilter *p_d = malloc(sizeof(DigitalFilter));
    if (!p_d) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    p_d->b_k = transform_polynomial(pa->order_numerator, -1.0, 1.0, 1.0, 1.0, pa->b_k);
    p_d->a_k = transform_polynomial(pa->order_denominator, -1.0, 1.0, 1.0, 1.0, pa->a_k);

    p_d->size_a = pa->size_a;
    p_d->size_b = pa->size_b;

    if (!p_d->b_k || !p_d->a_k) {
        free(p_d);
        fprintf(stderr, "Polynomial transformation failed.\n");
        return NULL;
    }

    return p_d;
}

void test_analog_to_digital()
{
    for (size_t order = 1; order < 10; order++) {
        for (size_t band = 0; band < 2; band++) {
            printf("Order %d Type %d Band %d \n", order, BUTTERWORTH, band);
            AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, band);

            if (!p) {
                fprintf(stderr, "Failed to generate the analog filter.\n");
                return;
            }

            DigitalFilter *p_d = transform_analog_to_digital(p);

            if (!p_d) {
                free_analog_filter(p);
                return;
            }

            normalize_to_b0(p_d);

            printf("\tNormalized Digital filter coefficients (a_k and b_k):\n");
            for (size_t i = 0; i < p_d->size_a; i++) {
                printf("\tak%zu %.32f\n", i, p_d->a_k[i]);
            }
            for (size_t i = 0; i < p_d->size_b; i++) {
                printf("\tbk%zu %.32f\n", i, p_d->b_k[i]);
            }

            free_analog_filter(p);
            free_digital_filter(p_d);
        }
    }
}
