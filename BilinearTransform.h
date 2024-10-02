#ifndef BILINEARTRANSFORM_H
#define BILINEARTRANSFORM_H

#include "AnalogFilters.h"
#include "DigitalFilters.h"

// Function prototypes for bilinear transformation and related operations

/**
 * Normalizes the filter coefficients to make the system causal.
 *
 * @param p Pointer to the DigitalFilter structure.
 * @param fs Sampling frequency.
 */
void horner_step6_make_causal_normalize_to_a0(DigitalFilter *p, long double fs);

/**
 * Shifts the polynomial coefficients with respect to a divisor.
 *
 * @param p Pointer to the DigitalFilter structure.
 * @param divisor Value to shift the polynomial.
 */
void horner_step2_5_shift_polynomial_with_n(DigitalFilter *p, long double divisor);

/**
 * Performs the first step of the Horner method by dividing and substituting.
 *
 * @param p Pointer to the DigitalFilter structure.
 * @param fs Sampling frequency.
 * @param prewarp_w0 Pre-warped frequency.
 * @param prewarp Boolean indicating whether to apply pre-warping.
 */
void horner_step1_divide_sn_substitute(DigitalFilter *p, long double fs, long double prewarp_w0, bool prewarp);

/**
 * Flips the coefficients of the filter.
 *
 * @param p Pointer to the DigitalFilter structure.
 */
void horner_step3_flip(DigitalFilter *p);

/**
 * Scales the polynomial zeros by a factor of 2.
 *
 * @param p Pointer to the DigitalFilter structure.
 */
void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p);

/**
 * The bilinear transformation compresses the entire left side of the s-plane
 * into a circle on the z-plane. This process is non-linear and emphasizes certain
 * parts of the frequency spectrum more than others.
 *
 * The emphasis depends on the sampling rate. At high sampling rates, the zeros
 * and poles corresponding to lower frequencies get pushed closer to z=1. This
 * behavior necessitates pre-warping.
 *
 * Low sampling rates can lead to aliasing and errors. Conversely, high sampling rates
 * require significant precision (more bits) to maintain accuracy, as critical information
 * concentrates near z=1. In this region, H(z) behaves normally, but high precision is needed
 * to differentiate values due to rounding and quantization errors.
 */

/**
 * Performs bilinear transformation using the Horner method.
 *
 * @param p Pointer to the AnalogFilter structure.
 * @param fs Sampling frequency.
 * @param w0 Pre-warped frequency.
 * @return Pointer to the resulting DigitalFilter structure.
 */
DigitalFilter *bilinear_transform_horner_method(AnalogFilter *p, double fs, double w0);

/**
 * Converts analog frequency to digital frequency.
 *
 * @param wa Analog frequency.
 * @param T Sampling period.
 * @return Digital frequency.
 */
double wa_2_wd(double wa, double T);

/**
 * Checks the stability of the system by examining reflection coefficients.
 *
 * @param A Coefficient array.
 * @param length Length of the coefficient array.
 * @return Stability status (1 for stable, 0 for unstable).
 */
int stabilitycheck(double *A, int length);

#endif // BILINEARTRANSFORM_H
