
#include "DigitalFilters.h"

void free_digital_filter(DigitalFilter *p)
{
    if (p != NULL) {
        free(p->a_k);
        free(p->b_k);
    }
}

void normalize_to_b0(DigitalFilter *p)
{
    double b_0 = p->a_k[0];
    for (size_t i = 0; i < p->size_a; i++) {
        p->a_k[i] *= (1.0 / b_0);
    }

    for (size_t i = 0; i < p->size_b; i++) {
        p->b_k[i] *= (1.0 / b_0);
    }
}
