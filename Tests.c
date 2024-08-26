
#include <stdlib.h>
#include <string.h>

#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "DigitalFilters.h"
#include "Filtering.h"
#include "Tests.h"§

const char *TEST_IMAGE_OUTPUT_PREFIX = "./tests/images/";

void test_analog_to_digital()
{
    for (size_t order = 1; order < 10; order++) {
        for (size_t band = 0; band < 2; band++) {
            printf("Order %zu Type %d Band %zu \n", order, BUTTERWORTH, band);
            AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, band);

            if (!p) {
                fprintf(stderr, "Failed to generate the analog filter.\n");
                return;
            }

            DigitalFilter *p_d = transform_analog_to_digital(p);

            if (!p_d) {
                free_analog_filter(p);
                return;
            }

            normalize_to_b0(p_d);

            printf("\tNormalized Digital filter coefficients (a_k and b_k):\n");
            for (size_t i = 0; i < p_d->size_a; i++) {
                printf("\tak%zu %.32f\n", i, p_d->a_k[i]);
            }
            for (size_t i = 0; i < p_d->size_b; i++) {
                printf("\tbk%zu %.32f\n", i, p_d->b_k[i]);
            }

            free_analog_filter(p);
            free_digital_filter(p_d);
        }
    }
}

void test_make_causal()
{
    for (size_t order = 1; order < 10; order++) {
        for (size_t band = 0; band < 2; band++) {
            printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

            AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, band);
            if (!p) {
                fprintf(stderr, "Failed to generate the analog filter.\n");
                return;
            }

            DigitalFilter *p_d = transform_analog_to_digital(p);
            if (!p_d) {
                free_analog_filter(p);
                return;
            }

            normalize_to_b0(p_d);

            DigitalFilterCausal *p_d_c = make_causal(p_d);
            if (!p_d_c) {
                free_analog_filter(p);
                free_digital_filter(p_d);
                return;
            }

            normalize_to_b0_causal(p_d_c);

            printf("\tNormalized Causal Digital filter coefficients (a_k and b_k):\n");
            for (size_t i = 0; i < p_d_c->size_a; i++) {
                printf("\tak%zu z^-%zu = %.32f\n", i, p_d_c->size_a - 1 - i, p_d_c->a_k[i]);
            }
            for (size_t i = 0; i < p_d_c->size_b; i++) {
                printf("\tak%zu z^-%zu = %.32f\n", i, p_d_c->size_b - 1 - i, p_d_c->b_k[i]);
            }

            free_analog_filter(p);
            free_digital_filter(p_d);
            free_causal_digital_filter(p_d_c);
        }
    }
}

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

void test_magnitude_phase_response_analog_digital()
{
    size_t n = 400;
    double *magnitudes_analog = (double *)calloc(n, sizeof(double));
    double *magnitudes_digital = (double *)calloc(n, sizeof(double));
    size_t order = 10;
    size_t band = 0;
    double cutoff = 0.5;

    printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

    AnalogFilter *p = generate_analog_filter(order, cutoff, BUTTERWORTH, band);
    if (!p) {
        fprintf(stderr, "Failed to generate the analog filter.\n");
        return;
    }

    DigitalFilter *p_d = transform_analog_to_digital(p);
    if (p_d == NULL) {
        free_analog_filter(p);
        return;
    }

    normalize_to_b0(p_d);

    DigitalFilterCausal *p_d_c = make_causal(p_d);
    if (p_d_c == NULL) {
        free_analog_filter(p);
        free_digital_filter(p_d);
        return;
    }

    normalize_to_b0_causal(p_d_c);

    char *l1 = concat_strings(4, "butter_analog_", "order_4_", "_lowpass", "color");
    char *l2 = concat_strings(4, "butter_digital_", "order_4_", "_lowpass", "color");
    char *labels[] = {l1,
                      l2};

    double *xs1 = magnitude_response_analog_filter(p, magnitudes_analog, n);
    double *xs2 = magnitude_response_causal_digital_filter(p_d_c, magnitudes_digital, n);

    char mag_filename[] = "analog_vs_digital.png";

    double *yss[] = {magnitudes_digital, magnitudes_digital};
    double *xss[] = {xs1,
                     xs1};

    if (plot_x_y_overlay(xss, yss, n, 2, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, mag_filename) != 0) {
        printf("png generation failed");
    }

    free(magnitudes_analog);
    free(magnitudes_digital);

    free(xs1);
    free(xs2);

    free_analog_filter(p);
    free_digital_filter(p_d);
    free_causal_digital_filter(p_d_c);
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
    plot_x_y_overlay(xss, yss, n, overlay_num, 4, labels, TEST_IMAGE_OUTPUT_PREFIX, filename);
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

void test_billinear_transform()
{
    size_t order = 10;
    AnalogFilter *p = generate_analog_filter(order, 1.0, 0, 0);
    // double coefficients_analog[] = {1.0000, 1.4142, 1.0000};
    for (size_t k = 0; k < p->size_a; k++)
        printf("an %f \n", p->a_k[k]);

    // double *transformed = transform_polynomial(order, -1.0, 1.0, 1.0, 1.0, p->a_k);
    double *r = (double *)calloc((order + 1), sizeof(double));
    memcpy(r, p->a_k, (order + 1) * sizeof(double));
    bilinear_transform(1.0, -1.0, 1.0, 1.0, p->a_k, order);

    for (size_t k = 0; k < order + 1; k++)
        printf("transform %20.17le \n", p->a_k[k]);

    free(r);
    free_analog_filter(p);
}



// Function to perform the bilinear transform
// implementation from the article is failed
// https://digitalcommons.unl.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=1085&context=imsefacpub
void bilinear_transform(double alpha, double beta, double delta, double gamma, double *r, int Order)
{
    int k, j;
    double sum;
    double *beta_powers, *gamma_powers;

    // Allocate memory for beta_powers and gamma_powers arrays
    beta_powers = (double *)malloc((Order + 1) * sizeof(double));
    gamma_powers = (double *)malloc((Order + 1) * sizeof(double));

    if (beta_powers == NULL || gamma_powers == NULL) {
        printf("Memory allocation failed!\n");
        return;
    }

    // Initialize the first elements of beta_powers and gamma_powers
    beta_powers[0] = gamma_powers[0] = 1;

    // Compute beta_powers and gamma_powers
    for (j = 1; j < Order; j++) {
        beta_powers[j] = beta * beta_powers[j - 1];
        gamma_powers[j] = gamma * gamma_powers[j - 1];
    }

    // Compute the coefficients using the bilinear transform
    for (k = 0; k < Order; k++) {
        sum = 0.0;
        for (j = 0; j < Order - k; j++) {
            sum += beta_powers[j] * gamma_powers[Order - k - j] * r[j];
        }
        for (j = 0; j < Order - k; j++) {
            r[j] = ((Order - k - j) * delta * r[j] + (j + 1) * alpha * r[j + 1]) / (k + 1);
        }
        r[Order - k] = sum;
    }

    // Free allocated memory
    free(beta_powers);
    free(gamma_powers);
}
