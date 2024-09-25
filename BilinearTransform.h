#ifndef BILINEARTRANSFORM_H
#define BILINEARTRANSFORM_H

#include "AnalogFilters.h"
#include "DigitalFilters.h"

void horner_step6_make_causal_normalize_to_a0(DigitalFilter *p, long double fs);

void horner_step2_5_shift_polynomial_with_n(DigitalFilter *p, long double divisor);

void horner_step1_divide_sn_substitute(DigitalFilter *p, long double fs, long double prewarp_w0, bool prewarp);

void horner_step3_flip(DigitalFilter *p);

void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p);

/*
A fast sampling rate isn't ideal for approximation.
The bilinear transformation compresses the entire left side of the s-plane into a circle on the z-plane.
This process is non-linear and emphasizes certain parts of the frequency spectrum more than others.

Which part is emphasized depends on the sampling rate.
At high sampling rates, the zeros and poles corresponding to lower frequencies get pushed closer to z=1.
This behavior is why pre-warping is necessary.

If the sample rate is low, aliasing will occur and cause errors. However, when the sample rate is high,
a significant amount of precision (more bits) is needed to meet accuracy requirements,
because much of the system's critical information will be concentrated near z=1. As a result,
H(z) will behave as if it's H(z)=1 almost everywhere except in regions close to z=1,
where poles and zeros are densely packed.

In that area, H(z) behaves more "normally," as designed, but high precision is needed to
differentiate between values in this narrow region. Essentially, the higher the sample rate,
the more challenging it becomes to maintain precision due to rounding and quantization errors.
*/

DigitalFilter *bilinear_transform_horner_method(AnalogFilter *p, double fs, double w0);

double wa_2_wd(double wa, double T);

int stabilitycheck(double *A, int length);

#endif
