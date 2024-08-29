#ifndef BILINEARTRANSFORM_H
#define BILINEARTRANSFORM_H

#include "AnalogFilters.h"
#include "DigitalFilters.h"

void horner_step3_flip(DigitalFilter *p);

void horner_step1_divide_sn_substitute(DigitalFilter *p, double fs);

void horner_shift_polynomial_with_n(DigitalFilter *p, double divisor);

void horner_step4_scale_polynomial_zeros_by_2(DigitalFilter *p);

DigitalFilter *bilinear_transform_horner_method(AnalogFilter *p, double fs);

void horner_step5_make_causal_normalize_to_b0(DigitalFilter *p);

#endif // BILINEARTRANSFORM_H
