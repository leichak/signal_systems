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
typedef struct
{
    double *y_n;
    double *x_n;
    double *b_k;
    double *a_k;
    size_t size_b;
    size_t size_a;
} FilterIIR;

FilterIIR *create_iir_filter(DigitalFilterCausal *p);

struct FilterFIR
{
    double *x_n;
    double *b_k;
};

#endif
