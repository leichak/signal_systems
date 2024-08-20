#ifndef BILINEARTRANSFORM_H
#define BILINEARTRANSFORM_H

#include "AnalogFilters.h"
#include "DigitalFilters.h"

/**
 * @brief Converts an analog filter to a digital filter using the bilinear transform technique.
 *
 * This function transforms an analog filter, specified by an `AnalogFilter` structure,
 * into an equivalent digital filter using the bilinear transformation method.
 *
 * @param p A pointer to the `AnalogFilter` structure representing the analog filter.
 * @return A pointer to the resulting `DigitalFilter` structure representing the digital filter.
 *
 * This function implements the bilinear transform method, which is commonly used to
 * convert continuous-time (analog) filters to discrete-time (digital) filters.
 *
 * Reference:
 *     For more information on the theory behind bilinear transform:
 *     https://digitalcommons.unl.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1085&context=imsefacpub
 */
DigitalFilter *transform_analog_to_digital(AnalogFilter *p);

/**
 * @brief Performs polynomial coefficient transformation as part of the bilinear transformation process.
 *
 * This function handles the transformation of polynomial coefficients during the bilinear
 * transform process, which is essential for converting analog filter coefficients to digital ones.
 *
 * @param Order The order of the polynomial to be transformed.
 * @param beta The beta coefficient used in the polynomial transformation.
 * @param gamma The gamma coefficient used in the polynomial transformation.
 * @param delta The delta coefficient used in the polynomial transformation.
 * @param alpha The alpha coefficient used in the polynomial transformation.
 * @param coefficients A pointer to the array containing the polynomial coefficients to be transformed.
 * @return A pointer to an array containing the transformed polynomial coefficients.
 *
 * The returned array is dynamically allocated and should be freed by the caller after it is no longer needed.
 */
double *transform_polynomial(int Order, double beta, double gamma, double delta, double alpha, double *coefficients);

#endif // BILINEARTRANSFORM_H
