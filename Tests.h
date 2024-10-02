#ifndef TESTS_H
#define TESTS_H

// Include necessary headers for the tests
#include "AnalogFilters.h"
#include "BilinearTransform.h"
#include "DigitalFilters.h"
#include "Filtering.h"
#include "FixedPoint.h"
#include "Utils.h"

/**
 * @brief Test function for validating the analog-to-digital filter transformation.
 *
 * This function tests the `transform_analog_to_digital()` function by generating an analog filter,
 * transforming it to a digital filter, and then outputting the filter coefficients.
 * It is used to ensure that the transformation process is working as expected.
 */
void test_analog_to_digital();

/**
 * @brief Function to plot frequency responses.
 *
 * This function visualizes the frequency responses of filters.
 */
void plot_freq_responses();

/**
 * @brief Test function for validating the DigitalFilter causal implementation.
 *
 * This function ensures that the causal implementation of digital filters behaves as expected.
 */
void test_make_causal();

/**
 * @brief Test function for generating various analog filters.
 *
 * This function is useful for validating filter generation logic.
 */
void test_generate_filters();

/**
 * @brief Test function for bilinear transformation.
 *
 * This function validates the correctness of the bilinear transform applied to filters.
 */
void test_bilinear_transform();

/**
 * @brief Test function for generating various filters of different orders.
 *
 * This function is useful for validating the logic behind filter generation for various orders.
 */
void test_various_orders_filters();

/**
 * @brief Test function for calculating magnitude and phase responses for analog and digital filters.
 *
 * This function is useful for validating the characteristics of filters.
 */
void test_magnitude_phase_response_analog_digital();

/**
 * @brief Test function for overlaying multiple lines on a plot.
 *
 * This function validates the correct overlaying of multiple data series.
 */
void test_overlay_multiple_lines();

/**
 * @brief Test function for generating color palettes.
 *
 * This function ensures that color generation for plots is functioning correctly.
 */
void test_generate_colors();

/**
 * @brief Test function for concatenating strings.
 *
 * This function validates the correctness of string concatenation operations.
 */
void test_concat_strings();

/**
 * @brief Test function for synthetic division.
 *
 * This function tests the implementation of synthetic division in polynomial arithmetic.
 */
void test_synthetic_division();

/**
 * @brief Test function for evaluating the error in a specified function.
 *
 * This function validates the correctness of the error evaluation in the EW function.
 */
void test_ew_function();

/**
 * @brief Test function for fixed-point multiplication.
 *
 * This function ensures that fixed-point multiplication operations are performed correctly.
 */
void test_fixed_multiplication();

/**
 * @brief Test function for quantization error with different k values.
 *
 * This function evaluates the quantization error associated with various fixed-point configurations.
 */
void test_quantization_error_different_k();

/**
 * @brief Test function for filtering with floating point numbers.
 *
 * This function validates the floating-point filtering process.
 */
void test_filtering_floating_point();

#endif // TESTS_H
