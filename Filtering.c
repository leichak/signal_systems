#include "Filtering.h"

DirectForm1 *create_df1(DigitalFilter *filter_ptr)
{
    DirectForm1 *df_ptr = (DirectForm1 *)malloc(sizeof(DirectForm1));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    df_ptr->y_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->x_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->max_size = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a : filter_ptr->size_b;

    for (size_t k = 0; k < df_ptr->size_a; k++)
        df_ptr->a_k[k] = filter_ptr->a_k[k];
    for (size_t k = 0; k < df_ptr->size_b; k++)
        df_ptr->b_k[k] = filter_ptr->b_k[k];
    return df_ptr;
}

DirectForm1Q15 *create_df1_q15(DigitalFilter *filter_ptr)
{
    DirectForm1Q15 *df_ptr = (DirectForm1Q15 *)malloc(sizeof(DirectForm1Q15));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (int *)calloc(filter_ptr->b_k, sizeof(int));
    df_ptr->a_k = (int *)calloc(filter_ptr->a_k, sizeof(int));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    df_ptr->y_k = (int *)calloc(filter_ptr->a_k, sizeof(int));
    df_ptr->x_k = (int *)calloc(filter_ptr->b_k, sizeof(int));
    df_ptr->max_size = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a : filter_ptr->size_b;

    for (size_t k = 0; k < df_ptr->size_a; k++)
        df_ptr->a_k[k] = (int)(to_fix(filter_ptr->a_k[k], 31));
    for (size_t k = 0; k < df_ptr->size_b; k++)
        df_ptr->b_k[k] = (int)(to_fix(filter_ptr->b_k[k], 31));

    return df_ptr;
}
DirectForm2 *create_df2(DigitalFilter *filter_ptr)
{
    DirectForm2 *df_ptr = (DirectForm2 *)malloc(sizeof(DirectForm2));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    size_t max_size = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a : filter_ptr->size_b;
    df_ptr->max_size = max_size;
    df_ptr->v_k = (double *)calloc(max_size, sizeof(double));

    for (size_t k = 0; k < df_ptr->size_a; k++)
        df_ptr->a_k[k] = filter_ptr->a_k[k];
    for (size_t k = 0; k < df_ptr->size_b; k++)
        df_ptr->b_k[k] = filter_ptr->b_k[k];

    return df_ptr;
}

DirectForm2Q15 *create_df2_q15(DigitalFilter *filter_ptr)
{
    DirectForm2Q15 *df_ptr = (DirectForm2Q15 *)malloc(sizeof(DirectForm2Q15));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    size_t max_size = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a : filter_ptr->size_b;
    df_ptr->max_size = max_size;
    df_ptr->v_k = (double *)calloc(max_size, sizeof(double));

    for (size_t k = 0; k < df_ptr->size_a; k++)
        df_ptr->a_k[k] = filter_ptr->a_k[k];
    for (size_t k = 0; k < df_ptr->size_b; k++)
        df_ptr->b_k[k] = filter_ptr->b_k[k];

    return df_ptr;
}

DirectForm1Transposed *create_df1_transposed(DigitalFilter *filter_ptr)
{
    DirectForm1Transposed *df_ptr = (DirectForm1Transposed *)malloc(sizeof(DirectForm1Transposed));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;

    size_t max_size = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a : filter_ptr->size_b;

    df_ptr->s_a = (double **)malloc(df_ptr->size_a * sizeof(double *));
    df_ptr->s_b = (double **)malloc(df_ptr->size_b * sizeof(double *));

    for (size_t k = 0; k < max_size; k++) { // each has 1 delay line so 2
        if (k < df_ptr->size_b)
            df_ptr->s_b[k] = (double *)calloc(2, sizeof(double));
        if (k < df_ptr->size_a)
            df_ptr->s_a[k] = (double *)calloc(2, sizeof(double));
    }

    return df_ptr;
}

DirectForm1TransposedQ15 *create_df1_transposed_q15(DigitalFilter *filter_ptr)
{
    DirectForm1TransposedQ15 *df_ptr = (DirectForm1TransposedQ15 *)malloc(sizeof(DirectForm1TransposedQ15));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (int *)calloc(filter_ptr->b_k, sizeof(int));
    df_ptr->a_k = (int *)calloc(filter_ptr->a_k, sizeof(int));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;

    df_ptr->s_a = (int *)calloc(2 * (df_ptr->size_b - 1), sizeof(int));
    df_ptr->s_a = (int *)calloc(2 * (df_ptr->size_a - 1), sizeof(int));

    return df_ptr;
}

DirectForm2Transposed *create_df2_transposed(DigitalFilter *filter_ptr)
{
    DirectForm2Transposed *df_ptr = (DirectForm2Transposed *)malloc(sizeof(DirectForm1Transposed));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    size_t max_order = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a - 1 : filter_ptr->size_b - 1;
    df_ptr->s_ab = (double *)calloc(max_order, sizeof(double));

    return df_ptr;
}

DirectForm2TransposedQ15 *create_df2_transposed_q15(DigitalFilter *filter_ptr)
{
    DirectForm2TransposedQ15 *df_ptr = (DirectForm2TransposedQ15 *)malloc(sizeof(DirectForm1TransposedQ15));

    if (df_ptr == NULL)
        return NULL;

    df_ptr->b_k = (double *)calloc(filter_ptr->b_k, sizeof(double));
    df_ptr->a_k = (double *)calloc(filter_ptr->a_k, sizeof(double));
    df_ptr->size_b = filter_ptr->size_b;
    df_ptr->size_a = filter_ptr->size_a;
    size_t max_order = filter_ptr->size_a > filter_ptr->size_b ? filter_ptr->size_a - 1 : filter_ptr->size_b - 1;
    df_ptr->s_ab = (double *)calloc(max_order, sizeof(double));

    return df_ptr;
}

void process_df1(DirectForm1 *ptr, double *x, double *y, size_t samples_num)
{
    for (size_t j = 0; j < samples_num; j++) {
        // Shifting x and new sample

        for (size_t k = 0; k < ptr->size_b - 1; k++)
            ptr->x_k[ptr->size_b - 1 - k] = ptr->x_k[ptr->size_b - 2 - k];
        ptr->x_k[0] = x[j];
        // Calculation of the newest
        double y_temp = 0.0;
        for (size_t k = 0; k < ptr->max_size; k++) {
            if (k < ptr->size_b)
                y_temp += ptr->x_k[k] * ptr->b_k[ptr->size_b - 1 - k];
            if (k > 0 && k < ptr->size_a)
                y_temp += ptr->y_k[k] * ptr->a_k[ptr->size_a - 1 - k];
        }
        // Shifting y and new sample
        for (size_t k = 0; k < ptr->size_a - 1; k++)
            ptr->y_k[ptr->size_a - 1 - k] = ptr->y_k[ptr->size_a - 2 - k];

        y[j] = (ptr->a_k[0] / ptr->b_k[0]) * y_temp;
        ptr->y_k[0] = y_temp;
    }
}

void process_df1_q15(DirectForm1Q15 *ptr, double *x, double *y, size_t samples_num)
{
    for (size_t j = 0; j < samples_num; j++) {
        // Shifting x and new sample
        for (size_t k = 0; k < ptr->size_b - 1; k++)
            ptr->x_k[ptr->size_b - 1 - k] = ptr->x_k[ptr->size_b - 2 - k];
        ptr->x_k[0] = to_fix(x[j], 15); // Q32,15 1 BIT FOR SIGN
        // Calculation of the newest
        int y_temp = 0;
        int mult_16_15 = 0;
        for (size_t k = 0; k < ptr->max_size; k++) {
            if (k < ptr->size_b) {
                mult_16_15 = fixed_N_k_mul(ptr->x_k[k], ptr->b_k[k], 15); // 2N, 2K, so shift
                y_temp = fixed_N_k_add(y_temp, mult_16_15);
            }
            if (k > 0 && k < ptr->size_a) {
                mult_16_15 = fixed_N_k_mul(ptr->y_k[k], ptr->a_k[k], 15); // 2N, 2K, so shift
                y_temp = fixed_N_k_add(y_temp, mult_16_15);
            }
        }
        // Shifting y and new sample
        for (size_t k = 0; k < ptr->size_a - 1; k++)
            ptr->y_k[ptr->size_a - 1 - k] = ptr->y_k[ptr->size_a - 2 - k];
        ptr->y_k[0] = y_temp;
        y[j] = (double)to_float(y_temp, 15);
    }
}

void process_df2(DirectForm2 *ptr, double *x, double *y, size_t samples_num)
{

    for (size_t j = 0; j < samples_num; j++) {
        // Calculation of the newest
        double v_temp = x[j];
        double y_temp = 0.0;
        for (size_t k = 0; k < ptr->max_size; k++) {
            if (k > 0 && k < ptr->size_a)
                v_temp += ptr->v_k[k] * ptr->a_k[k];
            if (k > 0 && k < ptr->size_b)
                y_temp += ptr->v_k[k] * ptr->b_k[k];
        }
        y_temp += v_temp * ptr->b_k[0];

        // Shifting y and new sample
        for (size_t k = 0; k < ptr->max_size - 1; k++)
            ptr->v_k[ptr->size_a - 1 - k] = ptr->v_k[ptr->size_a - 2 - k];
        ptr->v_k[0] = v_temp;
        y[j] = y_temp;
    }
}

void process_df2_q15(DirectForm2Q15 *ptr, double *x, double *y, size_t samples_num)
{

    for (size_t j = 0; j < samples_num; j++) {
        // Calculation of the newest
        int v_temp = to_fix(x[j], 15);
        int y_temp = 0;
        int mult_16_15 = 0;
        for (size_t k = 0; k < ptr->max_size; k++) {
            if (k > 0 && k < ptr->size_a) {
                mult_16_15 = fixed_N_k_mul(ptr->v_k[k], ptr->a_k[k], 15);
                v_temp = fixed_N_k_add(v_temp, mult_16_15);
            }
            if (k > 0 && k < ptr->size_b) {
                mult_16_15 = fixed_N_k_mul(ptr->v_k[k], ptr->b_k[k], 15);
                y_temp = fixed_N_k_add(y_temp, mult_16_15);
            }
        }
        mult_16_15 = fixed_N_k_mul(v_temp, ptr->b_k[0], 15);
        y_temp = fixed_N_k_add(y_temp, mult_16_15);

        // Shifting y and new sample
        for (size_t k = 0; k < ptr->max_size - 1; k++)
            ptr->v_k[ptr->size_a - 1 - k] = ptr->v_k[ptr->size_a - 2 - k];
        ptr->v_k[0] = v_temp;
        y[j] = y_temp;
    }
}

void process_df1_transposed(DirectForm1Transposed *ptr, double *x, double *y, size_t samples_num)
{
    // for (size_t j = 0; j < samples_num; j++) {

    //     // for (size_t k = 0; k < ptr->size_b - 1; k++)
    //     //     ptr->_[ptr->size_b - 1 - k] = ptr->x_k[ptr->size_b - 2 - k];
    //     // ptr->x_k[0] = x[j];
    //     double x_temp = x[j];

    //     for (size_t k = 1; k < ptr->size_a; k++) {
    //     }

    //     double y_temp = 0.0;
    //     for (size_t k = 0; k < ptr->max_size; k++) {
    //         if (k < ptr->size_b)
    //             y_temp += ptr->x_k[k] * ptr->b_k[k];
    //         if (k > 0 && k < ptr->size_a)
    //             y_temp += ptr->y_k[k] * ptr->a_k[k];
    //     }

    //     for (size_t k = 0; k < ptr->size_a - 1; k++)
    //         ptr->y_k[ptr->size_a - 1 - k] = ptr->y_k[ptr->size_a - 2 - k];
    //     ptr->y_k[0] = y[j];
    // }
}
void process_df1_transposed_q15(DirectForm1TransposedQ15 *ptr, double *x, double *y, size_t samples_num) {}
void process_df2_transposed(DirectForm2Transposed *ptr, double *x, double *y, size_t samples_num) {}
void process_df2_transposed_q15(DirectForm2TransposedQ15 *ptr, double *x, double *y, size_t samples_num) {}
