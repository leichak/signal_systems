#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AnalogFilters.h"
#include "DigitalFilters.h"

/**
 * @brief Performs the polynomial transformation required for the bilinear
 * transform.
 *
 * @param Order The order of the polynomial.
 * @param beta The beta coefficient used in the transformation.
 * @param gamma The gamma coefficient used in the transformation.
 * @param delta The delta coefficient used in the transformation.
 * @param alpha The alpha coefficient used in the transformation.
 * @param coefficients A pointer to the array of polynomial coefficients.
 * @return A pointer to an array containing the transformed polynomial
 * coefficients.
 */
double *transform_polynomial(int Order, double beta, double gamma, double delta,
                             double alpha, double *coefficients) {
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

  // Initialize the first powers to 1
  beta_powers[0] = gamma_powers[0] = 1.0;

  // Compute powers of beta and gamma
  for (int j = 1; j <= Order; j++) {
    beta_powers[j] = beta * beta_powers[j - 1];
    gamma_powers[j] = gamma * gamma_powers[j - 1];
  }

  // Perform the polynomial transformation
  for (int k = 0; k < Order; k++) {
    double sum = 0.0;

    for (int j = 0; j <= Order - k; j++) {
      coefficients[j] = ((Order - k - j) * delta * coefficients[j] +
                         (j + 1) * alpha * coefficients[j + 1]) /
                        (k + 1);
      sum += beta_powers[j] * gamma_powers[Order - k - j] * coefficients[j];
    }

    transformed_coefficients[Order - k] = sum;
  }

  // Process the last coefficient
  transformed_coefficients[0] = coefficients[0];

  // Free dynamically allocated memory for powers
  free(beta_powers);
  free(gamma_powers);

  // Return the array of transformed coefficients
  return transformed_coefficients;
}

/**
 * @brief Transforms an analog filter to a digital filter using bilinear
 * transformation.
 *
 * @param pa A pointer to the AnalogFilter struct.
 * @return A pointer to the DigitalFilter struct.
 */
DigitalFilter *transform_analog_to_digital(AnalogFilter *pa) {
  DigitalFilter *p_d = malloc(sizeof(DigitalFilter));

  if (!p_d) {
    fprintf(stderr, "Memory allocation failed.\n");
    return NULL;
  }

  p_d->b_k =
      transform_polynomial(pa->order_numerator, -1.0, 1.0, 1.0, 1.0, pa->b_k);
  p_d->a_k =
      transform_polynomial(pa->order_denominator, -1.0, 1.0, 1.0, 1.0, pa->a_k);
  p_d->size_a = pa->size_a;
  p_d->size_b = pa->size_b;

  if (p_d->b_k == NULL || p_d->a_k == NULL) {
    free(p_d);
    fprintf(stderr, "Polynomial transformation failed.\n");
    return NULL;
  }

  return p_d;
}

/**
 * @brief Test function for analog-to-digital filter transformation.
 */
void test_analog_to_digital() {
  int order = 5;
  AnalogFilter *p = generate_analog_filter(order, BUTTERWORTH);

  // Print analog coefficients
  printf("Analog filter coefficients (a_k) and (b_k):\n");
  for (size_t i = 0; i < p->size_a; i++) {
    printf("analog: ak%zu %.20f\n", i, p->a_k[i]);
  }
  for (size_t i = 0; i < p->size_b; i++) {
    printf("analog: bk%zu %.20f\n", i, p->b_k[i]);
  }

  if (!p) {
    fprintf(stderr, "Failed to generate analog filters.\n");
    free_analog_filter(p);
    return;
  }

  DigitalFilter *p_d = transform_analog_to_digital(p);

  if (!p_d) {
    free_digital_filter(p_d);
    return;
  }

  // Print digital coefficients
  printf("Digital filter coefficients (a_k and b_k):\n");
  for (size_t i = 0; i < p_d->size_a; i++) {
    printf("digital: ak%zu %.20f\n", i, p_d->a_k[i]);
  }
  for (size_t i = 0; i < p_d->size_b; i++) {
    printf("digital: bk%zu %.20f\n", i, p_d->b_k[i]);
  }

  normalize_to_b0(p_d);

  // Print digital coefficients
  printf("Digital filter coefficients (a_k and b_k):\n");
  for (size_t i = 0; i < p_d->size_a; i++) {
    printf("digital: ak%zu %.32f\n", i, p_d->a_k[i]);
  }
  for (size_t i = 0; i < p_d->size_b; i++) {
    printf("digital: bk%zu %.32f\n", i, p_d->b_k[i]);
  }

  // Free allocated memory
  free_analog_filter(p);
  free_digital_filter(p_d);
}