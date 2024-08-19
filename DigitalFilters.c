#include "DigitalFilters.h"

/**
 * @brief Frees the memory allocated for the digital filter structure.
 *
 * This function deallocates the memory used by the `DigitalFilter` structure,
 * including the memory allocated for the coefficients of the numerator
 * and denominator polynomials.
 *
 * @param p Pointer to the `DigitalFilter` structure to be freed.
 */
void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k); // Free memory allocated for a_k coefficients
        free(p->b_k); // Free memory allocated for b_k coefficients

        // Optional: nullify pointers to avoid dangling references
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

/**
 * @brief Normalizes the digital filter coefficients based on the first denominator coefficient (a0).
 *
 * This function normalizes the coefficients of the digital filter by dividing
 * all coefficients by the value of `a0` (the first coefficient of the
 * denominator polynomial). This is necessary to ensure the filter's scalability
 * with respect to `a0`.
 *
 * @param p Pointer to the `DigitalFilter` structure to be normalized.
 */
void normalize_to_b0(DigitalFilter *p)
{
    if (p == NULL)
        return; // Early return if the pointer is NULL

    double a_0 = p->a_k[0]; // Retrieve the first coefficient of the denominator polynomial

    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] /= a_0; // Normalize each denominator coefficient by a_0
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] /= a_0; // Normalize each numerator coefficient by a_0
    }
}
