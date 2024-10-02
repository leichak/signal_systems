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

// Generate colors based on the standard color palette
void generate_colors(int total_colors, RGBA *colors)
{
    for (int i = 0; i < total_colors; i++) {
        colors[i] = standard_colors[i % (sizeof(standard_colors) / sizeof(RGBA))];
    }
}

// Return the maximum of two integers
int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

// Return the maximum of two size_t values
size_t max_size_t(size_t a, size_t b)
{
    return (a > b) ? a : b;
}

// Fill an array with a range of values from start to stop
void fill_n_with_step(double *xs, size_t size, float start, float stop)
{
    float step = (stop - start) / (float)(size);
    xs[0] = start;
    for (size_t i = 1; i < size; i++) {
        xs[i] = xs[i - 1] + step;
    }
}

// Concatenate prefix and filename to create a full path
void *concat_path(char *path, char *prefix, char *filename)
{
    size_t prefix_len = strlen(prefix);
    size_t fname_len = strlen(filename);
    path = (char *)malloc((prefix_len + fname_len + 1) * sizeof(char)); // +1 for null terminator
    if (path == NULL)
        return NULL;

    strcpy(path, prefix);
    strcat(path, filename);
    return path;
}

// Create a scatter plot series with specified properties
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

// Plot multiple overlayed scatter plots
int plot_x_y_overlay(double *xs[], double *ys[], size_t n_points, size_t overlayed_num, size_t thickness, char *labels[], const char *prefix, char *filename, int range_y[])
{
    char *path = concat_path(NULL, prefix, filename);
    if (path == NULL)
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
    settings->autoPadding = true;

    wchar_t *xLabel = char_to_wchar(labels[overlayed_num]);
    settings->xLabel = xLabel;
    settings->xLabelLength = wcslen(xLabel);

    wchar_t *yLabel = char_to_wchar(labels[overlayed_num + 1]);
    settings->yLabel = yLabel;
    settings->yLabelLength = wcslen(yLabel);

    RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
    errorMessage = (StringReference *)malloc(sizeof(StringReference));

    success = DrawScatterPlotFromSettings(canvasReference, settings, errorMessage);
    if (!success) {
        fprintf(stderr, "Error: ");
        for (int i = 0; i < errorMessage->stringLength; i++) {
            fprintf(stderr, "%c", errorMessage->string[i]);
        }
        fprintf(stderr, "\n");
        free(path);
        return -1;
    }

    for (size_t k = 0; k < overlayed_num; k++) {
        wchar_t *label = char_to_wchar(labels[k]);
        DrawText(canvasReference->image, 10, 20 * (double)k, label, wcslen(label), GetBlack());
        free(label);
    }

    ByteArray *pngdata = ConvertToPNG(canvasReference->image);
    WriteToFile(pngdata, path);
    DeleteImage(canvasReference->image);

    FreeAllocations();
    free(path);

    return 0;
}

// Clip values in an array to a specified range
void clip_double(double *xs, size_t len, double max, double min)
{
    for (size_t k = 0; k < len; k++) {
        xs[k] = fmin(fmax(xs[k], min), max);
    }
}

// Concatenate multiple strings
char *concat_strings(int count, ...)
{
    va_list ap;
    va_start(ap, count);

    // Initialize with an empty string
    char *result = (char *)malloc(1);
    if (result == NULL) {
        va_end(ap);
        return NULL;
    }
    result[0] = '\0';

    for (int i = 0; i < count; i++) {
        char *next = va_arg(ap, char *);
        size_t needed = snprintf(NULL, 0, "%s %s", result, next) + 1; // +1 for null terminator
        result = realloc(result, needed);
        if (result == NULL) {
            va_end(ap);
            return NULL;
        }
        sprintf(result, "%s %s", result, next);
    }
    va_end(ap);

    return result;
}

// Convert a multibyte string to a wide character string
wchar_t *char_to_wchar(const char *str)
{
    if (str == NULL) {
        return NULL;
    }

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

// Unit conversion functions
double rad_s_2_hz(double omega)
{
    return omega * TWICE_PI_INVERSE;
}

double hz_2_rad_s(double f)
{
    return f * TWICE_PI;
}

double rad_s_2_rad_sample(double omega, double fs)
{
    return omega / fs;
}

double rad_sample_2_rad_s(double omega, double fs)
{
    return omega * fs;
}

double hz_2_cycle_sample(double f, double fs)
{
    return f / fs;
}

double cycle_sample_2_hz(double f, double fs)
{
    return f * fs;
}

double hz_2_half_cycle_sample(double f, double fs)
{
    return f / (fs / 2);
}

double half_cycle_sample_2_hz(double f, double fs)
{
    return f * (fs / 2);
}

// Calculate L2 norm mean
double l2_norm_mean(double *v1, double *v2, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return sqrt(sum) / size;
}

// Calculate L1 norm mean
double l1_norm_mean(double *v1, double *v2, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += fabs(v1[i] - v2[i]);
    }
    return sum / size;
}

// Calculate mean square error
double mean_square(double *v1, double *v2, int size)
{
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = v1[i] - v2[i];
        sum += diff * diff;
    }
    return sum / size;
}

// Generate a random float in a specified range
float random_float(float min, float max)
{
    float rand_value = (float)rand() / (float)RAND_MAX;
    return min + rand_value * (max - min);
}

// Generate a composite sine wave
double *generate_n_sines(double freqs[], size_t num, size_t samples_num, double fs)
{
    double *x = (double *)calloc(samples_num, sizeof(double));
    double A = 1.0 / (double)num;
    double pi_2 = TWICE_PI;

    for (size_t j = 0; j < samples_num; j++) {
        for (size_t k = 0; k < num; k++) {
            x[j] += A * sin((freqs[k] * (double)j / fs) * pi_2);
        }
    }

    return x;
}
