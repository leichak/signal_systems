#include "AnalogFilters.h"

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
    if (p != NULL)
    {
        free(p->a_k);
        free(p->b_k);
        free(p);
    }
}

void butter_sum_form_poles_coefficients(double *a_k, int length)
{
    a_k[0] = 1.0;
    double fi = M_PI / (2.0 * (length - 1));

    for (int i = 1; i < length; i++)
    {
        a_k[i] = a_k[i - 1] * cos((i - 1) * fi) / sin(i * fi);
    }
}

AnalogFilter *generate_analog_filter(int n)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    f->g_0 = 1.0;
    f->size_b = 1;
    f->b_k = (double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL)
    {
        free(f);
        return NULL;
    }
    f->b_k[0] = 1.0;

    f->size_a = n + 1;
    f->a_k = (double *)malloc(f->size_a * sizeof(double));
    if (f->a_k == NULL)
    {
        free(f->b_k);
        free(f);
        return NULL;
    }
    butter_sum_form_poles_coefficients(f->a_k, f->size_a);

    return f;
}

void transform_to_wc(AnalogFilter *p, double wc)
{
    for (int i = 1; i < p->size_a; i++)
    {
        p->a_k[i] *= pow(1.0 / wc, i);
    }
}

void normalize_to_max(AnalogFilter *p)
{
    double a_max = p->a_k[0];
    double b_max = p->b_k[0];

    for (int i = 1; i < p->size_b; i++)
    {
        if (p->b_k[i] > b_max)
            b_max = p->b_k[i];
    }

    for (int i = 1; i < p->size_a; i++)
    {
        if (p->a_k[i] > a_max)
            a_max = p->a_k[i];
    }

    p->g_0 = a_max / b_max;

    for (int i = 1; i < p->size_a; i++)
    {
        p->a_k[i] /= p->g_0;
    }
}

AnalogFilter *generate_analog_filter_wc(int n, double wc)
{
    AnalogFilter *f = (AnalogFilter *)malloc(sizeof(AnalogFilter));
    if (f == NULL)
        return NULL;

    f->g_0 = 1.0;
    f->size_b = 1;
    f->b_k = (double *)malloc(f->size_b * sizeof(double));
    if (f->b_k == NULL)
    {
        free(f);
        return NULL;
    }
    f->b_k[0] = 1.0;

    f->size_a = n + 1;
    f->a_k = (double *)malloc(f->size_a * sizeof(double));
    if (f->a_k == NULL)
    {
        free(f->b_k);
        free(f);
        return NULL;
    }
    butter_sum_form_poles_coefficients(f->a_k, f->size_a);
    transform_to_wc(f, wc);

    return f;
}

void test_coefficients_butter_sum_form_poles_coefficients()
{
    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order);
    AnalogFilter *p_wc = generate_analog_filter_wc(order, 0.2);

    for (size_t i = 0; i < p->size_a; i++)
    {
        printf("ak%zu %.4f ref: %.4f\n", i, p->a_k[i], ref5[i]);
    }

    free_analog_filter(p);
    free_analog_filter(p_wc);
}

void test_transform_wc()
{
    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order);
    transform_to_wc(p, 0.8);
    normalize_to_max(p);

    for (size_t i = 0; i < p->size_a; i++)
    {
        printf("g0 %.10f ak%zu %.10f ref: %.10f\n", p->g_0, i, p->a_k[i], ref5[i]);
    }

    free_analog_filter(p);
}
