#ifndef BILINEARTRANSFORM_H
#define BILINEARTRANSFORM_H

#include "AnalogFilters.h"
#include "DigitalFilters.h"

/**
 * @brief Performs a bilinear transform to convert an analog filter
 *        to a digital filter.
 *
 * This function transforms an analog filter specified by the
 * AnalogFilter structure to an equivalent digital filter using the
 * bilinear transform method.
 *
 * @param p A pointer to the AnalogFilter struct representing the analog filter.
 * @return A pointer to the resulting DigitalFilter struct representing the
 * digital filter.
 *
 * Reference:
 * https://digitalcommons.unl.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1085&context=imsefacpub
 */
DigitalFilter *transform_analog_to_digital(AnalogFilter *p);

/**
 * @brief Performs a polynomial transformation.
 *
 * This function handles the transformation of polynomial coefficients
 * as part of the bilinear transformation process.
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
                             double alpha, double *coefficients);

/**
 * @brief A test function to validate the analog to digital filter
 * transformation.
 *
 * This function is used to verify that the transform_analog_to_digital function
 * works correctly.
 */
void test_analog_to_digital();

#endif  // BILINEARTRANSFORM_H