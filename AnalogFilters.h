#ifndef ANALOGFILTERS_H
#define ANALOGFILTERS_H

#include "Utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Enumeration for filter types
typedef enum
{
    BUTTERWORTH = 0 // Currently, only Butterworth filter type is supported
} FilterTypes;

// Enumeration for various band types
typedef enum
{
    LOWPASS = 0,  // Lowpass filter type
    HIGHPASS = 1, // Highpass filter type
    BANDPASS = 2, // Bandpass filter type (not yet implemented)
    BANDSTOP = 3  // Bandstop (Notch) filter type (not yet implemented)
} BandType;

/**
 * @brief Calculate the linear gain of an nth-order Butterworth lowpass filter.
 *
 * @param n The order of the filter (number of poles).
 * @param g0 DC gain at zero frequency.
 * @param wc Cutoff frequency in radians per second (approximately -3dB).
 * @param w Frequency in radians per second at which the gain is calculated.
 * @return The linear gain of the filter.
 */
float butter_gain_low_pass(int n, float g0, float wc, float w);

/**
 * @brief Calculate the log gain of an nth-order Butterworth lowpass filter.
 *
 * @param n The order of the filter.
 * @param g0 DC gain at zero frequency (linear).
 * @param wc Cutoff frequency in radians per second (approximately -3dB).
 * @param w Frequency in radians per second at which the gain is calculated.
 * @return The log gain of the filter in decibels (dB).
 */
float butter_gain_low_pass_log10(int n, float g0, float wc, float w);

/**
 * @brief Structure to represent an analog filter.
 *
 * This structure is meant for representing analog filters (currently supports Butterworth filters).
 */
typedef struct
{
    long double *a_k;        ///< Coefficients of the denominator polynomial.
    long double *b_k;        ///< Coefficients of the numerator polynomial.
    double g_0;              ///< DC gain (gain at zero frequency).
    int size_a;              ///< Number of coefficients in the denominator.
    int size_b;              ///< Number of coefficients in the numerator.
    int order_numerator;     ///< Order of the numerator polynomial.
    int order_denominator;   ///< Order of the denominator polynomial.
    BandType band_type;      ///< Type of band for the filter.
    FilterTypes filter_type; ///< Type of filter (e.g., BUTTERWORTH).
} AnalogFilter;

/**
 * @brief Generates a Butterworth filter and initializes its coefficients.
 *
 * @param n The order of the filter (number of poles).
 * @param f Pointer to the AnalogFilter where the coefficients will be stored.
 * @return Pointer to the generated AnalogFilter structure.
 */
void *generate_butterworth(int n, AnalogFilter *f);

/**
 * @brief Free the memory allocated for an AnalogFilter structure.
 *
 * @param p Pointer to the AnalogFilter structure to be freed.
 */
void free_analog_filter(AnalogFilter *p);

/**
 * @brief Computes the Butterworth polynomial coefficients.
 *
 * This function calculates the coefficients of the Butterworth polynomial
 * using recursive and product formulas.
 *
 * @param order The order of the Butterworth filter.
 * @param coefficients Output array where the computed coefficients will be stored.
 *                     The size of the array should be `order + 1`.
 */
void butterworth_coefficients(int order, long double *coefficients);

/**
 * @brief Generate an analog filter prototype based on the specified filter type and band type.
 *
 * @param n Order of the filter (number of poles).
 * @param wc Cutoff frequency in radians per second.
 * @param filter_type Type of the filter (e.g., BUTTERWORTH).
 * @param band_type Type of the band (e.g., LOWPASS, HIGHPASS).
 * @return Pointer to an initialized AnalogFilter structure.
 */
AnalogFilter *generate_analog_filter(int n, double wc, FilterTypes filter_type, BandType band_type);

/**
 * @brief Transform the given filter into a low-pass filter with the specified cutoff frequency.
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param wc Desired cutoff frequency in radians per second.
 */
void transform_to_low_pass(AnalogFilter *p, double wc);

/**
 * @brief Transform the given low-pass filter into a high-pass filter with the specified cutoff frequency.
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param wc Desired cutoff frequency in radians per second.
 */
void transform_to_high_pass(AnalogFilter *p, double wc);

/**
 * @brief Normalize the coefficients of the numerator and denominator by their maximum value.
 *
 * @param p Pointer to the AnalogFilter structure.
 */
void normalize_to_max(AnalogFilter *p);

/**
 * @brief Calculate the magnitude response of a causal filter.
 *
 * M(\omega)=|H(e^{j\omega})| is a function of \omega called the magnitude response of the amplitude response of the system.
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param magnitudes Output array where the magnitude response will be stored.
 * @param n The number of frequency points.
 * @return Array of frequencies (w_k).
 */
double *magnitude_response_analog_filter(AnalogFilter *p, double *magnitudes, int n);

/**
 * @brief Calculate the phase response of a causal filter.
 *
 * Phase response is computed by solving the equation:
 * \tan(\theta(\omega)) = \frac{\Im\{H(e^{j\omega})\}}{\Re\{H(e^{j\omega})\}}
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param phases Output array where the phase response will be stored.
 * @param n The number of frequency points.
 * @return Array of frequencies (w_k).
 */
double *phase_response_analog_filter(AnalogFilter *p, double *phases, int n);

#endif // ANALOGFILTERS_H
