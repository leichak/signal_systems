#include "AnalogFilters.h"

#include <assert.h>

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

AnalogFilter *generate_analog_filter(int n, FilterTypes filter_type)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    switch (filter_type) {
    case BUTTERWORTH:
        assert(generate_butterworth(n, f));
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
    butterworth_coefficients(n, f->a_k);
}

void transform_to_wc(AnalogFilter *p, double wc)
{
    for (int i = 1; i < p->size_a; i++) {
        p->a_k[i] *= pow(1.0 / wc, i);
    }
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

AnalogFilter *generate_analog_filter_wc(int n, double wc,
                                        FilterTypes filter_type)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    f = generate_analog_filter(n, filter_type);
    transform_to_wc(f, wc);

    return f;
}

void test_coefficients_butter_sum_form_poles_coefficients()
{
    int order = 10;
    AnalogFilter *p = generate_analog_filter(order, BUTTERWORTH);

    if (p == NULL) {
        printf("Failed to allocate");
    }

    for (size_t i = 0; i < p->size_a; i++) {
        printf("ak%zu %.4f ref: %.4f\n", i, p->a_k[i], p->a_k[i]);
    }
    free_analog_filter(p);
}

void test_transform_wc()
{
    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order, BUTTERWORTH);
    transform_to_wc(p, 0.8);
    normalize_to_max(p);

    for (size_t i = 0; i < p->size_a; i++) {
        printf("g0 %.10f ak%zu %.10f ref: %.10f\n", p->g_0, i, p->a_k[i], ref5[i]);
    }

    free_analog_filter(p);
}
