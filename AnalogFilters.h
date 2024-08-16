#ifndef ANALOGFILTERS_H
#define ANALOGFILTERS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Calculate the linear gain of an nth-order Butterworth lowpass filter.
 *
 * @param n The order of the filter (number of poles).
 * @param g0 DC gain at zero frequency.
 * @param wc Cutoff frequency (approximately -3dB).
 * @param w Frequency in radians per second.
 * @return The linear gain of the filter.
 */
float butter_gain_low_pass(int n, float g0, float wc, float w);

/**
 * @brief Calculate the log gain of an nth-order Butterworth lowpass filter.
 *
 * @param n The order of the filter.
 * @param g0 DC gain at zero frequency (linear).
 * @param wc Cutoff frequency (approximately -3dB).
 * @param w Frequency in radians per second.
 * @return The log gain of the filter in decibels (dB).
 */
float butter_gain_low_pass_log10(int n, float g0, float wc, float w);

/**
 * @brief Structure to represent an analog Butterworth filter.
 */
typedef struct
{
    double *a_k; ///< Coefficients of the denominator polynomial.
    double *b_k; ///< Coefficients of the numerator polynomial.
    double g_0;  ///< DC gain.
    int size_a;  ///< Number of coefficients in the denominator.
    int size_b;  ///< Number of coefficients in the numerator.
} AnalogFilter;

/**
 * @brief Free the memory allocated for an AnalogFilter structure.
 *
 * @param p Pointer to the AnalogFilter structure to be freed.
 */
void free_analog_filter(AnalogFilter *p);

/**
 * @brief Recursively generate the coefficients for the Butterworth polynomial in sum form.
 *
 * @param a_k Pointer to the array where coefficients will be stored.
 * @param length Number of coefficients to generate.
 */
void butter_sum_form_poles_coefficients(double *a_k, int length);

/**
 * @brief Generate an analog Butterworth filter prototype.
 *
 * @param n Order of the filter (number of poles).
 * @return Pointer to an AnalogFilter structure.
 */
AnalogFilter *generate_analog_filter(int n);

/**
 * @brief Transform the Butterworth filter to have a specified cutoff frequency.
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param wc Desired cutoff frequency.
 */
void transform_to_wc(AnalogFilter *p, double wc);

/**
 * @brief Normalize the coefficients of the numerator and denominator to a maximum value.
 *
 * @param p Pointer to the AnalogFilter structure.
 */
void normalize_to_max(AnalogFilter *p);

/**
 * @brief Generate an analog Butterworth filter prototype with a specified cutoff frequency.
 *
 * @param n Order of the filter (number of poles).
 * @param wc Cutoff frequency in radians per second.
 * @return Pointer to an AnalogFilter structure.
 */
AnalogFilter *generate_analog_filter_wc(int n, double wc);

/**
 * @brief Test generated coefficients for Butterworth analog filter.
 */
void test_coefficients_butter_sum_form_poles_coefficients();

/**
 * @brief Test the transformation of Butterworth coefficients with a specified cutoff frequency.
 */
void test_transform_wc();

#endif // ANALOGFILTERS_H
