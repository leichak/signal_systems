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

/// @brief Represent analog butter filter
typedef struct
{
    double *a;
    int size;
} ButterworthAnalogNormalizedPrototypePolynomialFactors;

/// @brief Generate prototype analog filter
/// @param n
/// @return
ButterworthAnalogNormalizedPrototypePolynomialFactors generate(int n)
{
    ButterworthAnalogPrototypePolynomialFactors prototype;

    return prototype;
}
