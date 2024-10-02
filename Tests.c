#include "Tests.h"
#include <complex.h>
#include <stdlib.h>
#include <string.h>

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

// Test magnitude and phase response for analog and digital filters
void test_magnitude_phase_response_analog_digital()
{
    size_t n = 2000;
    double *magnitudes_analog = (double *)calloc(n, sizeof(double));
    double *magnitudes_digital = (double *)calloc(n, sizeof(double));

    size_t order = 2;
    size_t band = LOWPASS;
    double cutoff = 0.1 * M_PI; // pi = fs/2, 2pi = fs
    double fs = 1.0;

    printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

    AnalogFilter *p = generate_analog_filter(order, cutoff, BUTTERWORTH, band);
    if (!p) {
        fprintf(stderr, "Failed to generate the analog filter.\n");
        return;
    }

    DigitalFilter *p_d = bilinear_transform_horner_method(p, fs, cutoff);
    if (p_d == NULL) {
        free_analog_filter(p);
        return;
    }

    for (size_t k = 0; k < p_d->size_a; k++)
        printf("a%zu: %.64f\n", k, p_d->a_k[k]);

    char *l1 = concat_strings(4, "butter_analog_", "order_4_", "_lowpass", "color");
    char *l2 = concat_strings(4, "butter_digital_", "order_4_", "_lowpass", "color");
    char *labels[] = {l1, l2, "Frequency - w", "Magnitude"};

    double *xs1 = magnitude_response_analog_filter(p, magnitudes_analog, n);
    double *xs2 = magnitude_response_digital_filter(p_d, magnitudes_digital, n, fs);

    char mag_filename[] = "analog_vs_digital.png";

    // Clipping
    clip_double(magnitudes_digital, n, 1.0, -96.0);

    double *yss[] = {magnitudes_analog, magnitudes_digital};
    double *xss[] = {xs1, xs1};

    int y_range[] = {-96.0, 1.0};
    if (plot_x_y_overlay(xss, yss, n, 2, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, mag_filename, y_range) != 0) {
        printf("PNG generation failed");
    }

    horner_step3_flip(p_d);
    printf("Stable: %d \n", stabilitycheck(p_d->a_k, p_d->size_a));

    // Free allocated memory
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
    double ys5[] = {-10.0, -0.3, 0.5, 0.0, -3.0, -1.0};
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
        .a_k = &a_k, // Coefficients of the denominator polynomial.
        .b_k = &b_k, // Coefficients of the numerator polynomial.
        .size_a = 3, // Number of coefficients in the denominator polynomial.
        .size_b = 1, // Number of coefficients in the numerator polynomial.
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

    double *a_k = (double *)calloc(3, sizeof(double));
    double *b_k = (double *)calloc(1, sizeof(double));

    for (size_t k = 0; k < 1; k++) {
        b_k[k] = b[k];
    }
    for (size_t k = 0; k < 3; k++) {
        a_k[k] = a[k];
    }

    DigitalFilter p_s = {
        .a_k = a_k,  // Coefficients of the denominator polynomial.
        .b_k = b_k,  // Coefficients of the numerator polynomial.
        .size_a = 3, // Number of coefficients in the denominator polynomial.
        .size_b = 1, // Number of coefficients in the numerator polynomial.
        .power_numerator = 0,
        .power_denominator = 0,
    };

    DigitalFilter *p = &p_s;
    double fs = 1.0;

    // Display coefficients before substitution
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }

    printf("\nSubstituting and dividing...\n");
    horner_step1_divide_sn_substitute(p, fs, 0.0, false);

    // Display coefficients after substitution
    for (size_t k = 0; k < p->size_a; k++) {
        printf("\t a%zu %Lf", k, p->a_k[k]);
    }
    for (size_t k = 0; k < p->size_b; k++) {
        printf("\t b%zu %Lf", k, p->b_k[k]);
    }

    free(a_k);
    free(b_k);
}

void test_free()
{
    DigitalFilter *p = (DigitalFilter *)malloc(sizeof(DigitalFilter));
    if (p) {
        free(p);
    }
}
