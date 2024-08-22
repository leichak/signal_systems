
#include <stdlib.h>
#include <string.h>

#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "DigitalFilters.h"
#include "Filtering.h"
#include "PbPlots/pbPlots.h"
#include "PbPlots/supportLib.h"
#include "Tests.h"

static char *TEST_IMAGE_OUTPUT_PREFIX = "./tests/images/";

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

void test_magnitude_response()
{
    size_t n = 100;
    double *mag_analog = (double *)calloc(n, sizeof(double));
    double *mag_digital = (double *)calloc(n, sizeof(double));
    size_t order = 4;
    size_t band = 0;

    printf("Order %zu Type %d Band %zu\n", order, BUTTERWORTH, band);

    AnalogFilter *p = generate_analog_filter(order, 0.5, BUTTERWORTH, band);
    if (!p) {
        fprintf(stderr, "Failed to generate the analog filter.\n");
        return;
    }

    // DigitalFilter *p_d = transform_analog_to_digital(p);
    // if (!p_d) {
    //     free_analog_filter(p);
    //     return;
    // }

    // normalize_to_b0(p_d);

    // DigitalFilterCausal *p_d_c = make_causal(p_d);
    // if (!p_d_c) {
    //     free_analog_filter(p);
    //     free_digital_filter(p_d);
    //     return;
    // }

    // normalize_to_b0_causal(p_d_c);

    // printf("\tNormalized Causal Digital filter coefficients (a_k and b_k):\n");
    // for (size_t i = 0; i < p_d_c->size_a; i++) {
    //     printf("\tak%zu z^-%zu = %.32f\n", i, p_d_c->size_a - 1 - i, p_d_c->a_k[i]);
    // }
    // for (size_t i = 0; i < p_d_c->size_b; i++) {
    //     printf("\tak%zu z^-%zu = %.32f\n", i, p_d_c->size_b - 1 - i, p_d_c->b_k[i]);
    // }

    double *xs1 = magnitude_response_analog_filter(p, mag_analog, n);
    // double *xs2 = magnitude_response_causal_digital_filter(p_d_c, mag_digital, n);

    // for (size_t k = 0; k < n; k++) {
    //     printf("%f \n", mag_analog[k]);
    // }

    if (plot_x_y_single(xs1, mag_analog, n) != 0) {
        printf("pdf generation failed");
    }

    free(mag_analog);
    free(mag_digital);
    free(xs1);
    // free(xs2);

    free_analog_filter(p);
    // free_digital_filter(p_d);
    // free_causal_digital_filter(p_d_c);
}

int plot_x_y_single(double *xs, double *ys, size_t n)
{
    char filename[] = "magnitudes.png";
    size_t prefix_len = strlen(TEST_IMAGE_OUTPUT_PREFIX);
    size_t fname_len = strlen(&filename);
    char path[prefix_len + fname_len];
    strcpy(path, TEST_IMAGE_OUTPUT_PREFIX);
    strcat(path, &filename);

    _Bool success;

    StartArenaAllocator();

    RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
    StringReference *errorMessage = CreateStringReference(L"", 0);

    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = xs;
    series->xsLength = n;
    series->ys = ys;
    series->ysLength = n;
    series->linearInterpolation = true;
    series->lineType = L"solid";
    series->lineTypeLength = wcslen(series->lineType);
    series->lineThickness = 2;
    series->color = GetGray(0.3);
    ScatterPlotSeries *s[] = {series};

    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 1920;
    settings->height = 1080;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->scatterPlotSeries = s;
    settings->scatterPlotSeriesLength = 1;

    success = DrawScatterPlotFromSettings(canvasReference, settings, errorMessage);

    if (success) {
        size_t length;
        ByteArray *pngdata = ConvertToPNG(canvasReference->image);
        WriteToFile(pngdata, path);
        DeleteImage(canvasReference->image);
    } else {
        fprintf(stderr, "Error: ");
        for (int i = 0; i < errorMessage->stringLength; i++) {
            fprintf(stderr, "%c", errorMessage->string[i]);
        }
        fprintf(stderr, "\n");
    }

    FreeAllocations();

    return success ? 0 : 1;
}
