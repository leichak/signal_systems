#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
Algorithm implement:
    It finds the low-pass analog prototype poles, zeros, and gain using the function buttap.
    It converts the poles, zeros, and gain into state-space form.
    If required, it uses a state-space transformation to convert the lowpass filter into a bandpass, highpass, or bandstop filter with the desired frequency constraints.
    For digital filter design, it uses bilinear to convert the analog filter into a digital filter through a bilinear transformation with frequency prewarping.
    Careful frequency adjustment enables the analog filters and the digital filters to have the same frequency response magnitude at Wn or at w1 and w2.
    It converts the state-space filter back to its transfer function or zero-pole-gain form, as required.
*/

/// @brief The linear gain of nth-order butterworth lowpass filter
/// @param n order of filter (linear)
/// @param g0 dc gain at zero frequency
/// @param wc cutoff frequency approx. -3dB
/// @param w frequency in rad/s
/// @return
float butter_gain_low_pass(int n, float g0, float wc, float w)
{
    return (pow(g0, 2) / (1 + pow((w / wc), 2 * n)));
}

/// @brief The log gain of nth-order butterworth lowpass filter
/// @param n order of filter
/// @param g0 dc gain at zero frequency (linear)
/// @param wc cutoff frequency approx. -3dB
/// @param w frequency in rad/s
/// @return
float butter_gain_low_pass_log10(int n, float g0, float wc, float w)
{
    return 10.0 * log10((pow(g0, 2) / (1 + pow((w / wc), 2 * n))));
}

/// @brief Represent analog filter in a sum form
/// derived from: https://en.wikipedia.org/wiki/Butterworth_filter#Normalized_Butterworth_polynomials
typedef struct
{
    double *a_k; // denominator
    double *b_k; // numerator
    int size_a;  // size of numerator
    int size_b;  // size of denominator
} AnalogFilter;

/// @brief freeing AnalogFilter struct
/// @param p struct pointer
void free_analog_filter(AnalogFilter *p)
{
    free(p->a_k);
    free(p->b_k);
    free(p);
}

/// @brief Function that recursively generates butter sum form a_k coefficient
/// accordingly to https://en.wikipedia.org/wiki/Butterworth_filter#Normalized_Butterworth_polynomials
/// @param a_k pointer to coefficients array
/// @param length length for an array
void butter_sum_form_poles_coefficients(double *a_k, int length)
{
    a_k[0] = 1.0;
    double fi = M_PI / (2.0 * (double)(length - 1));
    size_t k;
    for (size_t i = 1; i < length; i++)
    {
        k = i;
        a_k[k] = a_k[k - 1] * cos((double)(k - 1) * fi) / sin((double)(k)*fi);
    }
}

/// @brief Generate prototype analog filter, for now butter
/// @param n order of a filter (number of poles)
/// @return pointer to a struct
AnalogFilter *generate_analog_filter(int n)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    // nominator of normalized is always size of 1 and is 1
    f->size_b = 1;
    f->b_k = (double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL)
        return NULL;
    f->b_k[0] = 1.0;
    f->size_a = n + 1;
    f->a_k = (double *)malloc(f->size_a * sizeof(double));
    if (f->a_k == NULL)
        return NULL;
    butter_sum_form_poles_coefficients(f->a_k, f->size_a);

    return f;
}

/// Tests

/// @brief Test generated coefficients for butter analog filter
void test_coefficients_buffer_form_poles_coefficients()
{

    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order);

    for (size_t i = 0; i < p->size_a; i++)
    {
        printf("ak%zu %.4f ref: %.4f \n", i, p->a_k[i], ref5[i]);
    }

    // create some asserts
}
