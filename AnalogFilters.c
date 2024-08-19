#include "AnalogFilters.h"

#include <assert.h>

int max_int(int a, int b)
{
    if (a > b) {
        return a;
    } else
        return b;
}

float butter_gain_low_pass(int n, float g0, float wc, float w)
{
    return (pow(g0, 2) / (1 + pow((w / wc), 2 * n)));
}

float butter_gain_low_pass_log10(int n, float g0, float wc, float w)
{
    return 10.0f * log10(pow(g0, 2) / (1 + pow((w / wc), 2 * n)));
}

void free_analog_filter(AnalogFilter *p)
{
    if (p != NULL) {
        // Free the dynamically allocated coefficients of the denominator polynomial
        // (a_k)
        if (p->a_k != NULL) {
            free(p->a_k);
            p->a_k = NULL; // Prevent accidental use after free
        }

        // Free the dynamically allocated coefficients of the numerator polynomial
        // (b_k)
        if (p->b_k != NULL) {
            free(p->b_k);
            p->b_k = NULL; // Prevent accidental use after free
        }

        // Free the AnalogFilter structure itself
        free(p);
        p = NULL; // Optional step, nullifying local pointer passed by value
    }
}

void butterworth_coefficients(int order, double *coefficients)
{
    // Initialize γ (gamma)
    double gamma = M_PI / (2.0 * order);

    // Initialize the first coefficient (a0)
    coefficients[0] = 1.0;

    // Compute remaining coefficients using the product formula
    for (int k = 1; k <= order; ++k) {
        // Use the recursive relation with the product formula to compute a_k
        coefficients[k] =
            coefficients[k - 1] * (cos((k - 1) * gamma) / sin(k * gamma));
    }

    // Due to symmetry a_k = a_n-k
    for (int k = 0; k <= order / 2; ++k) {
        coefficients[order - k] = coefficients[k];
    }
}

AnalogFilter *generate_analog_filter(int n, double wc, FilterTypes filter_type, BandType band_type)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    switch (filter_type) {
    case BUTTERWORTH:
        assert(generate_butterworth(n, f) != NULL);
    }

    switch (band_type) {
    case LOWPASS:
        transform_to_low_pass(f, wc);
        break;
    case HIGHPASS:
        transform_to_high_pass(f, wc);
        break;
    default:
        transform_to_low_pass(f, wc);
    }

    return f;
}

void *generate_butterworth(int n, AnalogFilter *f)
{
    f->g_0 = 1.0;
    f->order_denominator = n;
    f->order_numerator = 0;
    f->size_b = f->order_numerator + 1;
    f->b_k = (double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL) {
        free(f);
        return NULL;
    }
    f->b_k[0] = 1.0;

    f->size_a = n + 1;
    f->a_k = (double *)malloc(f->size_a * sizeof(double));
    if (f->a_k == NULL) {
        free(f->b_k);
        free(f);
        return NULL;
    }
    butterworth_coefficients(f->order_denominator, f->a_k);

    return f; // meaning it worked
}

void transform_to_low_pass(AnalogFilter *p, double wc)
{
    // Substitute s -> s / wc
    for (int i = 1; i < p->size_a; i++) {
        p->a_k[i] *= pow(1.0 / wc, i);
    }
}

void transform_to_high_pass(AnalogFilter *p, double wc)
{
    // Substitute s -> wc / s,
    // Multiplying by highest orders

    // Saving orders
    int order_numerator = p->order_numerator;
    int order_denominator = p->order_denominator;

    printf("%d %d \n", order_denominator, order_numerator);

    // Transformation && wc
    for (size_t i = 1; i < p->size_a; i++) {
        p->a_k[i] *= pow(wc, i);
    }
    for (size_t i = 1; i < p->size_b; i++) {
        p->b_k[i] *= pow(wc, i);
    }

    // Now changing order by numerator/denominator order multiplication (shifting)
    p->order_denominator += order_numerator;
    p->size_a = p->order_denominator + 1;
    p->order_numerator += order_denominator;
    p->size_b = p->order_numerator + 1;
    double *b = (double *)malloc(p->size_b * sizeof(double));
    double *a = (double *)malloc(p->size_a * sizeof(double));
    size_t size = max_int(p->size_a, p->size_b);

    for (size_t i = 0; i < size; i++) {
        if (i < p->size_b) {
            b[i] = 0.0;
            if (i >= order_denominator)
                b[i] = p->b_k[i - order_denominator];
        }
        if (i < p->size_a) {
            a[i] = 0.0;
            if (i >= order_numerator) {
                a[i] = p->a_k[i - order_numerator];
            }
        }
    }

    // Reversing orders (multiplying by polynomial order)
    int i = 0, j = p->size_a - 1;
    while (i < j || i != j) {
        double t = a[j];
        a[j] = a[i];
        a[i] = t;
        i++, j--;
    }
    i = 0, j = p->size_b - 1;
    while (i < j || i != j) {
        double t = b[j];
        b[j] = b[i];
        b[i] = t;
        i++, j--;
    }

    free(p->a_k);
    free(p->b_k);
    p->a_k = a;
    p->b_k = b;
}

void normalize_to_max(AnalogFilter *p)
{
    double a_max = p->a_k[0];
    double b_max = p->b_k[0];

    for (int i = 1; i < p->size_b; i++) {
        if (p->b_k[i] > b_max)
            b_max = p->b_k[i];
    }

    for (int i = 1; i < p->size_a; i++) {
        if (p->a_k[i] > a_max)
            a_max = p->a_k[i];
    }

    p->g_0 = a_max / b_max;

    for (int i = 1; i < p->size_a; i++) {
        p->a_k[i] /= p->g_0;
    }
}

void test_generate_filters()
{
    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, HIGHPASS);

    for (size_t i = 0; i < p->size_a; i++) {
        printf("g0 %.10f ak%zu %.20f ref: %.20f\n", p->g_0, i, p->a_k[i], ref5[i]);
    }

    free_analog_filter(p);
}
