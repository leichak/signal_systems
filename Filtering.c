#include "Filtering.h"

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
