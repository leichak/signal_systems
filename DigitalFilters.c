#include "DigitalFilters.h"
#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "Utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

void free_causal_digital_filter(DigitalFilterCausal *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
        p->a_k = NULL;
        p->b_k = NULL;
    }
}

void normalize_to_b0(DigitalFilter *p)
{
    if (p == NULL) {
        return;
    }

    double a_0 = p->a_k[0];

    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] /= a_0;
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] /= a_0;
    }
}

void normalize_to_b0_causal(DigitalFilterCausal *p)
{
    if (p == NULL) {
        return;
    }

    double a_0 = p->a_k[p->size_a - 1];

    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] /= a_0;
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] /= a_0;
    }
}

DigitalFilterCausal *make_causal(DigitalFilter *p)
{
    if (!p) {
        return NULL;
    }

    DigitalFilterCausal *p_d_c = (DigitalFilterCausal *)malloc(sizeof(DigitalFilterCausal));
    assert(p_d_c != NULL);

    size_t max_size = p->size_a > p->size_b ? p->size_a : p->size_b;
    p_d_c->a_k = (double *)calloc(max_size, sizeof(double));
    p_d_c->b_k = (double *)calloc(max_size, sizeof(double));
    assert(p_d_c->a_k != NULL && p_d_c->b_k != NULL);

    memcpy(p_d_c->a_k, p->a_k, sizeof(double) * p->size_a);
    memcpy(p_d_c->b_k, p->b_k, sizeof(double) * p->size_b);

    p_d_c->size_a = max_size;
    p_d_c->size_b = max_size;

    return p_d_c;
}

void frequency_response_causal_digital_filter(DigitalFilterCausal *p, double *freq, int n)
{
}
