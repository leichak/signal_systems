
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "Tests.h"

const char *TEST_IMAGE_OUTPUT_PREFIX = "./tests/images/";

// Test function for generating analog filters
void test_generate_filters()
{
    int order = 4;
    double ref5[5] = {1.0, 2.6131, 3.4142, 2.6131, 1.0};
    AnalogFilter *p = generate_analog_filter(order, 1.0, BUTTERWORTH, LOWPASS);

    for (size_t i = 0; i < p->size_a; i++) {
        printf("g0 %.10f a_k%zu %.20f ref: %.20f\n", p->g_0, i, p->a_k[i], ref5[i]);
    }

    free_analog_filter(p);
}

void test_various_orders_filters()
{
    for (size_t order = 1; order < 10; order++) {
        for (size_t band = 0; band < 2; band++) {
            printf("Order %zu Type %d Band %zu \n", order, BUTTERWORTH, band);
            AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, band);
            for (size_t i = 0; i < p->size_a; i++) {
                printf("\tg0 %.10f a_k%zu %.20f \n", p->g_0, i, p->a_k[i]);
            }
            free_analog_filter(p);
        }
    }
}
// What are inputs to both plots analog and digital - normalized freq
// Plot stability of filters after Billinear transform
// It works for particularly low frequencies
// So stability might check whether poles are inside unit circle
// because for higher orders higher fs cause weird distortions
void test_magnitude_phase_response_analog_digital()
{
    size_t n = 2000;
    double *magnitudes_analog = (double *)calloc(n, sizeof(double));
    double *magnitudes_digital = (double *)calloc(n, sizeof(double));
    size_t order = 1;
    size_t band = LOWPASS;
    double cutoff = 0.2 * M_PI; // pi = fs/2, 2pi = fs
    double fs = 48000.0;

    printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

    AnalogFilter *p = generate_analog_filter(order, cutoff, BUTTERWORTH, band);
    if (!p) {
        fprintf(stderr, "Failed to generate the analog filter.\n");
        return;
    }

    DigitalFilter *p_d = bilinear_transform_horner_method(p, fs, cutoff); // why
    if (p_d == NULL) {
        free_analog_filter(p);
        return;
    }

    for (size_t k = 0; k < p_d->size_a; k++)
        printf("a%zu: %f\n", k, p_d->a_k[k]);

    char *l1 = concat_strings(4, "butter_analog_", "order_4_", "_lowpass", "color");
    char *l2 = concat_strings(4, "butter_digital_", "order_4_ ", "_lowpass", "color");
    char *labels[] = {l1,
                      l2, "Frequency - w", "Magnitude"};

    double *xs1 = magnitude_response_analog_filter(p, magnitudes_analog, n);
    double *xs2 = magnitude_response_digital_filter(p_d, magnitudes_digital, n, fs);

    char mag_filename[] = "analog_vs_digital.png";

    // Clipping
    clip_double(magnitudes_digital, n, 1.0, -96.0);

    double *yss[] = {magnitudes_analog, magnitudes_digital};
    double *xss[] = {xs1,
                     xs1};

    int y_range[] = {-96.0, 1.0};
    if (plot_x_y_overlay(xss, yss, n, 2, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, mag_filename, y_range) != 0) {
        printf("png generation failed");
    }

    free(magnitudes_analog);
    free(magnitudes_digital);

    free(xs1);
    free(xs2);

    free_analog_filter(p);
    free_digital_filter(p_d);
}

void test_overlay_multiple_lines()
{
    char filename[] = "overlayed1.png";
    char *labels[] = {(char *)""};

    size_t n = 6;
    double ys1[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double ys2[] = {2.0, 1.0, 0.5, 0.0, 1.0, -1.0};
    double ys3[] = {5.0, 1.0, 1.5, 0.0, -2.0, -1.0};
    double ys4[] = {5.0, -0.0, -1.5, 0.0, -2.0, 1.0};
    double ys5[] = {-10.0, -.3, 0.5, 0.0, -3.0, -1.0};
    double xs[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double *xss[] = {xs, xs, xs, xs, xs};
    double *yss[] = {ys1, ys2, ys3, ys4, ys5};

    size_t overlay_num = 5;
    int y_range[] = {1, -50};

    plot_x_y_overlay(xss, yss, n, overlay_num, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, filename, y_range);
}

void test_generate_colors()
{
    RGBA colors[10];
    generate_colors(10, colors);
}

void test_concat_strings()
{
    char *concatenated = concat_strings(4, "label1", "label2", "label3", "label4");

    printf("%s", concatenated);
}

void test_synthetic_division()
{
    double a_k[] = {1.0, 2.0, 2.0}; // powers 0,1,2,3....
    double b_k[] = {2.0};

    DigitalFilter p = {
        .a_k = &a_k, ///< Coefficients of the denominator polynomial.
        .b_k = &b_k, ///< Coefficients of the numerator polynomial.
        .size_a = 3, ///< Number of coefficients in the denominator polynomial.
        .size_b = 1, ///< Number of coefficients in the numerator polynomial.
        .power_numerator = 0,
        .power_denominator = 2,
    };
    horner_step2_5_shift_polynomial_with_n(&p, 1.0);

    for (size_t k = 0; k < p.size_a; k++) {
        printf("\t a%zu %f", k, p.a_k[k]);
    }
    for (size_t k = 0; k < p.size_b; k++) {
        printf("\t b%zu %f", k, p.b_k[k]);
    }

    printf("\n");
}

void test_bilinear_transform()
{
    double a[] = {2.0, 2.0, 1.0};
    double b[] = {2.0};
    double *a_k = (double *)calloc(3, sizeof(double)); //{2.0, 2.0, 1.0}; // powers 0,1,2,3....
    double *b_k = (double *)calloc(1, sizeof(double));
    for (size_t k = 0; k < 1; k++)
        b_k[k] = b[k];
    for (size_t k = 0; k < 3; k++)
        a_k[k] = a[k];

    DigitalFilter p_s = {
        .a_k = a_k,  ///< Coefficients of the denominator polynomial.
        .b_k = b_k,  ///< Coefficients of the numerator polynomial.
        .size_a = 3, ///< Number of coefficients in the denominator polynomial.
        .size_b = 1, ///< Number of coefficients in the numerator polynomial.c
        .power_numerator = 0,
        .power_denominator = 0,
    };

    DigitalFilter *p = &p_s;
    double fs = 1.0;

    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\nsubstitute / divide\n");
    horner_step1_divide_sn_substitute(p, fs, 0.0, false);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\n");
    printf("\tRef: N (s) → N (x) = 0.5x2 D(s) → D(x) = 0.5x2 + x + 1 \n");
    printf("\nshift + 1\n");
    horner_step2_5_shift_polynomial_with_n(p, 1);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\n");
    printf("\tRef: N (x + 1) = 0.5x2 + x + 0.5 D(x + 1) = 0.5x2 + 2x + 2.5 \n");
    printf("\nflip\n");
    horner_step3_flip(p);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\n");
    printf("\tRef: N ( 1/x + 1) = 0.5x2 + x + 0.5 D( 1/x + 1) = 2.5x2 + 2x + 0.5 \n");
    printf("\nscale by 2\n ");
    horner_step4_scale_polynomial_zeros_by_2(p);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\n");
    printf("\tRef: N ( 2/x + 1) = 0.125x2 + 0.5x + 0.5 D( 2/x + 1) = 0.625x2 + x + 0.5 \n");
    printf("\nshift -1\n");
    horner_step2_5_shift_polynomial_with_n(p, -1);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }
    printf("\n");
    printf("\tRef: N (x) = 0.125x2 + 0.25x + 0.125 D(x) = 0.625x2 − 0.25x + 0.125 \n");
    printf("\nmake causal and normalise\n");
    horner_step6_make_causal_normalize_to_a0(p, 1);
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t ak%d %Lf", p->power_denominator + k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t bk%d %Lf", p->power_numerator + k, p->b_k[k]);
    }
    printf("\n\tRef: Y [z](1 − 0.4z−1 + 0.2z−2) = X[z](0.2 + 0.4z−1 + 0.2z−2)\n");
}

void test_ew_function()
{
    char filename[] = "ew_test.png";
    char *labels[] = {(char *)"real", (char *)"imag", (char *)"abs"};

    size_t n = 10;
    size_t overlay_num = 3;

    double fs = 20;
    double T = 1 / fs;
    double *x0 = (double *)calloc(n, sizeof(double));
    double *y0real = (double *)calloc(n, sizeof(double));
    double *y0imag = (double *)calloc(n, sizeof(double));
    double *y0abs = (double *)calloc(n, sizeof(double));
    double *y0mag = (double *)calloc(n, sizeof(double));
    fill_n_with_step(x0, n, -M_PI, M_PI);
    for (size_t k = 0; k < n; k++) {
        y0real[k] = creal(cexp(2 * M_PI * k / (double)n * I));
        y0imag[k] = cimag(cexp(2 * M_PI * k / (double)n * I));
        y0abs[k] = cabs(cexp(M_PI * k / (double)n * I));
    }
    double *xss[] = {x0, y0real, x0};
    double *yss[] = {y0real, y0imag, y0abs};

    int y_range[] = {1, -50};

    plot_x_y_overlay(xss, yss, n, overlay_num, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, filename, y_range);

    free(x0);
    for (size_t p = 0; p < 3; p++) {
        free(yss[p]);
    }
}

void test_fixed_multiplication()
{
    float c[5] = {0.2, -0.4, 0.8, -0.4, 0.2};
    float taps[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    float result;

    int c_8_7[5];
    int taps_8_7[5];
    int result_16_14 = 0;

    size_t i;

    // convert c float to c_8_7 (fix<8,7>), taps to fix,
    for (i = 0; i < 5; i++) {
        c_8_7[i] = (int)(c[i] * (1 << 7)); // multiply 2^k
        taps_8_7[i] = (int)(taps[i] * (1 << 7));
        result_16_14 += c_8_7[i] * taps_8_7[i];
    }

    // Output data type is fix<16,14>
    printf("Multiplication  result %f \n", (result_16_14 * 1.0f / (1 << 14)));
}

void test_quantization_error_different_k()
{
    size_t size = 1000;
    double random_numbers[size];
    for (size_t k = 0; k < size; k++)
        random_numbers[k] = (double)random_float(-1.0, 1.0);

    // Different k
    size_t ks_num = 15;
    double ks[ks_num];
    double errors[ks_num];
    for (size_t k = 1; k <= ks_num; k++)
        ks[k - 1] = k;

    // Fill results vectors with random numbers
    double **results = malloc(sizeof(double *) * ks_num);
    for (size_t k = 1; k <= ks_num; k++) {
        results[k - 1] = (double *)malloc(size * sizeof(double));
        memcpy(results[k - 1], &random_numbers, size * sizeof(double));
        for (size_t j = 0; j < size; j++) // convert and convert back
            results[k - 1][j] = to_float(to_fix(results[k - 1][j], k), k);
        errors[k - 1] = 20.0 * log10(l1_norm_mean(results[k - 1], &random_numbers, size) / 1.0);
    }

    char filename[] = "k_fractional_vs_error.png";
    char *labels[] = {(char *)"error", (char *)"number of fractional bits", (char *)"20*log10(l1_norm_error)"};

    size_t n = ks_num;
    size_t overlay_num = 1;

    double **xss[] = {&ks};
    double **yss[] = {&errors};

    int range[] = {1, -50};

    plot_x_y_overlay(xss, yss, n, overlay_num, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, filename, range);

    for (size_t k = 0; k < ks_num; k++)
        free(results[k]);
}

void test_filtering_floating_point()
{

    size_t order = 3;
    size_t band = LOWPASS;
    double cutoff = 0.3 * M_PI / 2.0; // 400 hz
    double fs = 1.0;

    printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

    AnalogFilter *p = generate_analog_filter(order, cutoff, BUTTERWORTH, band);
    if (!p) {
        fprintf(stderr, "Failed to generate the analog filter.\n");
        return;
    }

    DigitalFilter *p_d = bilinear_transform_horner_method(p, fs, cutoff); // why
    if (p_d == NULL) {
        free_analog_filter(p);
        return;
    }

    double A[] = {0.1, 0.1, -0.4};
    int stable = stabilitycheck(p_d->a_k, p_d->size_a);

    printf("Stable %d \n", stable);

    DirectForm1 *p_df_1 = create_df1(p_d);
    if (p_df_1 == NULL) {
        free_analog_filter(p);
        free_digital_filter(p_d);
        return;
    }

    // Generate sines
    double freqs[] = {30.0, 900.0};
    size_t parts_num = 1;
    double time_len = 0.5;
    size_t samples_num = (size_t)(fs * time_len);
    double *x = generate_n_sines(freqs, parts_num, fs, time_len);
    double *y = (double *)calloc(samples_num, sizeof(double));
    double *time = (double *)calloc(samples_num, sizeof(double));
    fill_n_with_step(time, samples_num, 0.0, time_len);

    // Processing
    process_df1(p_df_1, x, y, samples_num);

    char filename[] = "df_floating_filtering.png";
    char *labels[] = {(char *)"original", (char *)"filtered", (char *)"time", (char *)"amp"};

    size_t overlay_num = 2;

    for (size_t k = 0; k < 1; k++)
        printf("%f \n", y[k]);

    double *xss[] = {time, time};
    double *yss[] = {x, y};

    int range[] = {1, -50};
    plot_x_y_overlay(xss, yss, 300, overlay_num, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, filename, range);

    free(x);
    free(time);
    free(y);
    free_analog_filter(p);
    free_digital_filter(p_d);
}
