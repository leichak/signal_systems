#ifndef FILTERING_H
#define FILTERING_H

#include "DigitalFilters.h"
#include "Utils.h"

/**
 * Assuming that filter is causal in the form of Z transform implement filtering
 * reference:
 * https://www.mathworks.com/help/signal/ug/filter-implementation-and-analysis.html
 *
 */
struct FilterIIR
{
    double *y_n;
    double *x_n;
    double *b_k;
    double *a_k;
    size_t size_b;
    size_t size_a;
};

FilterIIR *create_iir_filter(DigitalFilterCausal *p)
{
    FilterIIR *p_f_iir = (FilterIIR *)malloc(sizeof(FilterIIR));

    if (p_f_iir == NULL) {
        return NULL;
    }

    p_f_iir->a_k = (double *)calloc(p->size_a - 1, sizeof(double));
    if (p_f_iir->a_k == NULL) {
        free(p_f_iir);
        return NULL;
    }
    p_f_iir->b_k = (double *)calloc(p->size_b, sizeof(double));
    if (p_f_iir->b_k == NULL) {
        free(p_f_iir);
        return NULL;
    }
    p_f_iir->y_n = (double *)calloc((p->size_a - 1), sizeof(double));
    if (p_f_iir->y_n == NULL) {
        free(p_f_iir);
        return NULL;
    }
    p_f_iir->x_n = (double *)calloc((p->size_a), sizeof(double));
    if (p_f_iir->x_n == NULL) {
        free(p_f_iir);
        return NULL;
    }
    p_f_iir->size_a = p->size_a - 1;
    p_f_iir->size_b = p->size_b;

    size_t size = max_int(p->size_a - 1, p->size_b);
    for (size_t i = 0; i < size; i++) {
        if (i < p_f_iir->size_b)
            p_f_iir->b_k[i] = p->b_k[i];
        if (i < p_f_iir->size_a)
            p_f_iir->a_k[i] = p->a_k[i];
    }

    return p_f_iir;
}

struct FilterFIR
{
    double *x_n;
    double *b_k;
};

#endif
