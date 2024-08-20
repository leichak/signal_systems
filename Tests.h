#ifndef TESTS_H
#define TESTS_H

// Include pb plots
#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "DigitalFilters.h"
#include "Filtering.h"

/**
 * @brief Test function for validating the analog-to-digital filter transformation.
 *
 * This function tests the `transform_analog_to_digital()` function by generating an analog filter,
 * transforming it to a digital filter, and then outputting the filter coefficients. It is used
 * to ensure that the transformation process is working as expected.
 */
void test_analog_to_digital();

void plot_freq_responses();

/**
 * @brief this function is testing DigitalFilter Causal implementation
 */
void test_make_causal();

/**
 * @brief Test function for generating various analog filters.
 *
 * This function is useful for validating filter generation logic.
 */
void test_generate_filters();

/**
 * @brief Test function for generating various analog filters.
 *
 * This function is useful for validating filter generation logic.
 */
void test_various_orders_filters();

#endif
