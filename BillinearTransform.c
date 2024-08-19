#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AnalogFilters.h"
#include "DigitalFilters.h"

/**
 * @brief Performs the polynomial transformation required for the bilinear transform.
 *
 * This transformation is typically used in converting an analog filter's transfer function
 * to its digital counterpart.
 *
 * @param Order The order of the polynomial.
 * @param beta The beta coefficient used in the transformation.
 * @param gamma The gamma coefficient used in the transformation.
 * @param delta The delta coefficient used in the transformation.
 * @param alpha The alpha coefficient used in the transformation.
 * @param coefficients A pointer to the array of polynomial coefficients to be transformed.
 * @return A pointer to an array containing the transformed polynomial coefficients.
 */
double *transform_polynomial(int Order, double beta, double gamma, double delta, double alpha, double *coefficients)
{
    // Allocate memory for the arrays used in the transformation
    double *beta_powers = malloc((Order + 1) * sizeof(double));
    double *gamma_powers = malloc((Order + 1) * sizeof(double));
    double *transformed_coefficients = malloc((Order + 1) * sizeof(double));

    // Check for memory allocation failures
    if (!beta_powers || !gamma_powers || !transformed_coefficients) {
        fprintf(stderr, "Memory allocation failed.\n");
        free(beta_powers);
        free(gamma_powers);
        free(transformed_coefficients);
        return NULL;
    }

    // Initialize the first powers of beta and gamma to 1
    beta_powers[0] = gamma_powers[0] = 1.0;

    // Compute the powers of beta and gamma for the transformation
    for (int j = 1; j <= Order; j++) {
        beta_powers[j] = beta * beta_powers[j - 1];
        gamma_powers[j] = gamma * gamma_powers[j - 1];
    }

    // Perform the polynomial transformation
    for (int k = 0; k < Order; k++) {
        double sum = 0.0;
        for (int j = 0; j <= Order - k; j++) {
            coefficients[j] = ((Order - k - j) * delta * coefficients[j] + (j + 1) * alpha * coefficients[j + 1]) / (k + 1);
            sum += beta_powers[j] * gamma_powers[Order - k - j] * coefficients[j];
        }
        transformed_coefficients[Order - k] = sum;
    }

    // Handle the last coefficient
    transformed_coefficients[0] = coefficients[0];

    // Free dynamically allocated memory
    free(beta_powers);
    free(gamma_powers);

    // Return the array containing transformed coefficients
    return transformed_coefficients;
}

/**
 * @brief Transforms an analog filter to a digital filter using bilinear transformation.
 *
 * This function converts the analog Butterworth filter into its digital equivalent.
 *
 * @param pa Pointer to the AnalogFilter structure.
 * @return Pointer to the DigitalFilter structure.
 */
DigitalFilter *transform_analog_to_digital(AnalogFilter *pa)
{
    // Allocate memory for the DigitalFilter structure
    DigitalFilter *p_d = malloc(sizeof(DigitalFilter));
    if (!p_d) {
        fprintf(stderr, "Memory allocation failed.\n");
        return NULL;
    }

    // Perform polynomial transformation for numerator and denominator
    p_d->b_k = transform_polynomial(pa->order_numerator, -1.0, 1.0, 1.0, 1.0, pa->b_k);
    p_d->a_k = transform_polynomial(pa->order_denominator, -1.0, 1.0, 1.0, 1.0, pa->a_k);

    // Store sizes of the polynomial coefficients
    p_d->size_a = pa->size_a;
    p_d->size_b = pa->size_b;

    // Check if the transformation has failed
    if (!p_d->b_k || !p_d->a_k) {
        free(p_d);
        fprintf(stderr, "Polynomial transformation failed.\n");
        return NULL;
    }

    return p_d;
}

/**
 * @brief Test function for analog-to-digital filter transformation.
 *
 * This function performs a test by generating an analog filter and transforming it
 * into a digital filter using the bilinear transformation method.
 */
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

            // Perform the analog to digital transformation
            DigitalFilter *p_d = transform_analog_to_digital(p);

            if (!p_d) {
                free_analog_filter(p);
                return;
            }
            // Normalize coefficients by b_0
            normalize_to_b0(p_d);

            // Print the normalized digital filter coefficients
            printf("\tNormalized Digital filter coefficients (a_k and b_k):\n");
            for (size_t i = 0; i < p_d->size_a; i++) {
                printf("\tak%zu %.32f\n", i, p_d->a_k[i]);
            }
            for (size_t i = 0; i < p_d->size_b; i++) {
                printf("\tbk%zu %.32f\n", i, p_d->b_k[i]);
            }
            // Free allocated memory
            free_analog_filter(p);
            free_digital_filter(p_d);
        }
    }
}
