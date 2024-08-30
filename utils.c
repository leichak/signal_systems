#include "Utils.h"

#include <locale.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <wchar.h>

static char *TEST_IMAGE_OUTPUT_PREFIX = "./tests/images/";

const RGBA standard_colors[] = {
    {1.0, 0.0, 0.0, 1.0},          // Red
    {0.0, 1.0, 0.0, 1.0},          // Green
    {0.0, 0.0, 1.0, 1.0},          // Blue
    {1.0, 1.0, 0.0, 1.0},          // Yellow
    {0.0, 1.0, 1.0, 1.0},          // Cyan
    {1.0, 0.0, 1.0, 1.0},          // Magenta
    {0.7529, 0.7529, 0.7529, 1.0}, // Silver (192/255)
    {0.5019, 0.5019, 0.5019, 1.0}, // Gray (128/255)
    {0.0, 0.0, 0.0, 1.0},          // Black
    {1.0, 1.0, 1.0, 1.0}           // White
};

void generate_colors(int total_colors, RGBA *colors)
{
    int i;
    for (i = 0; i < total_colors; i++) {
        colors[i] = standard_colors[i % total_colors];
    }
}

int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

size_t max_size_t(size_t a, size_t b)
{
    return (a > b) ? a : b;
}

void fill_n_with_step(double *xs, size_t size, float start, float stop)
{
    float step = (stop - start) / (float)(size);
    xs[0] = start;
    for (int i = 1; i < size; i++) {
        xs[i] = xs[i - 1] + step;
    }
}

void *concat_path(char *path, char *prefix, char *filename)
{
    size_t prefix_len = strlen(prefix);
    size_t fname_len = strlen(filename);
    path = (char *)malloc((prefix_len + fname_len) * sizeof(char));
    if (path == NULL)
        return NULL;
    strcpy(path, prefix);
    strcat(path, filename);
    return path;
}

ScatterPlotSeries *get_series_solid_thickness_color(double *xs, double *ys, size_t n, size_t thickness, RGBA *color)
{
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
    series->xs = xs;
    series->xsLength = n;
    series->ys = ys;
    series->ysLength = n;
    series->linearInterpolation = true;
    series->lineType = L"solid";
    series->lineTypeLength = wcslen(series->lineType);
    series->lineThickness = thickness;
    series->color = color;

    return series;
}

int plot_x_y_overlay(double *xs[], double *ys[], size_t n_points, size_t overlayed_num, size_t thickness, char *labels[], const char *prefix, char *filename)
{

    char *path;
    path = concat_path(path, prefix, filename);

    if (concat_path(path, prefix, filename) == NULL)
        return 1;

    _Bool success;

    StartArenaAllocator();

    StringReference *errorMessage;

    ScatterPlotSeries *s[overlayed_num];
    RGBA colors[overlayed_num];
    generate_colors(overlayed_num, colors);

    for (size_t k = 0; k < overlayed_num; k++) {
        s[k] = get_series_solid_thickness_color(xs[k], ys[k], n_points, thickness, &colors[k]);
    }

    ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
    settings->width = 1920;
    settings->height = 1080;
    settings->scatterPlotSeries = s;
    settings->scatterPlotSeriesLength = overlayed_num;
    settings->autoBoundaries = true;
    settings->autoPadding = false;

    RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
    errorMessage = (StringReference *)malloc(sizeof(StringReference));

    success = DrawScatterPlotFromSettings(canvasReference, settings, errorMessage);
    if (!success) {
        return -1;
    }

    for (size_t k = 0; k < overlayed_num; k++) {
        wchar_t *label = char_to_wchar(labels[k]);
        DrawText(canvasReference->image, 10, 20 * (double)k, label, wcslen(label), GetBlack());
        free(label);
        if (!success)
            return -1;
    }

    if (success) {
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

    free(path);

    return success ? 0 : 1;
}

char *concat_strings(int count, ...)
{
    va_list ap;
    va_start(ap, count);
    char *prev = "";
    for (int i = 0; i < count; i++) {
        char *next = va_arg(ap, char *);
        size_t needed = snprintf(NULL, 0, "%s %s", prev, next);
        char *buffer = malloc(needed);
        sprintf(buffer, "%s %s", prev, next);
        prev = buffer;
    }
    va_end(ap);

    return prev;
}

wchar_t *char_to_wchar(const char *str)
{
    if (str == NULL) {
        return NULL;
    }

    // Set the locale to the user's default to handle multibyte to wide character conversion
    setlocale(LC_CTYPE, "");

    // Calculate the length of the wide string
    size_t len = mbstowcs(NULL, str, 0);
    if (len == (size_t)-1) {
        perror("mbstowcs");
        return NULL;
    }

    // Allocate memory for the wide string
    wchar_t *wstr = malloc((len + 1) * sizeof(wchar_t));
    if (wstr == NULL) {
        perror("malloc");
        return NULL;
    }

    // Convert the multibyte string to a wide string
    mbstowcs(wstr, str, len + 1);

    return wstr;
}

double rad_s_2_hz(double omega)
{
    return omega * TWICE_PI_INVERSE;
}

double rad_s_2_hz_inv(double f)
{
    return f * TWICE_PI;
}

double rad_s_2_rad_sample(double omega, double fs)
{
    return omega / fs;
}

double rad_s_2_rad_sample_inv(double omega, double fs)
{
    return omega * fs;
}

double hz_2_cycle_sample(double f, double fs)
{
    return f / fs;
}

double hz_2_cycle_sample_inv(double f, double fs)
{
    return f * fs;
}

double hz_2_half_cycle_sample(double f, double fs)
{
    return f / (fs / 2);
}

double hz_2_half_cycle_sample_inv(double f, double fs)
{
    return f * (fs / 2);
}
