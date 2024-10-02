#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "PbPlots/pbPlots.h"
#include "PbPlots/supportLib.h"

// Constants
#define TWICE_PI (2.0 * M_PI)
#define TWICE_PI_INVERSE (1.0 / (2.0 * M_PI))

// String Utilities
void *concat_path(char *path, char *prefix, char *filename);
char *concat_strings(int count, ...);
wchar_t *char_to_wchar(const char *str);

// Drawing Utilities
void HSLToRGB(double h, double s, double l, unsigned char *r, unsigned char *g, unsigned char *b);
ScatterPlotSeries *get_series_solid_thickness_color(double *xs, double *ys, size_t n, size_t thickness, RGBA *color);
int plot_x_y_single(double *xs, double *ys, size_t n, const char *prefix, char *filename);
int plot_x_y_overlay(double *xs[], double *ys[], size_t n_points, size_t overlayed_num, size_t thickness, char *labels[], const char *prefix, char *filename, int range[]);
void generate_colors(int m, RGBA *colors);

// Pre-computed Utilities
void fill_n_with_step(double *xs, size_t size, float start, float stop);
void clip_double(double *xs, size_t len, double max, double min);

// Numerical Utilities
int max_int(int a, int b);

// Units Utilities

/// @brief Convert rad/s to Hz. 1 Hz equals 2π rad/s.
/// @param omega Omega in PI radians/s.
/// @return Frequency in Hz (cycles/s).
double rad_s_2_hz(double omega);

/// @brief Convert Hz to rad/s. 1 Hz equals 2π rad/s.
/// @param f Frequency in Hz (cycles/s).
/// @return Omega in PI radians/s.
double hz_2_rad_s(double f);

/// @brief Convert rad/s to rad/sample.
/// @param omega Omega in PI radians/s.
/// @param fs Sampling frequency in samples/s.
/// @return Radians/sample.
double rad_s_2_rad_sample(double omega, double fs);

/// @brief Convert rad/sample to rad/s.
/// @param omega Radians/sample.
/// @param fs Sampling frequency in samples/s.
/// @return Omega in PI radians/s.
double rad_sample_2_rad_s(double omega, double fs);

/// @brief Convert Hz to cycles/sample.
/// @param f Frequency in cycles/second.
/// @param fs Sampling frequency in samples/second.
/// @return Cycles/sample.
double hz_2_cycle_sample(double f, double fs);

/// @brief Convert cycles/sample to Hz.
/// @param f Frequency in cycles/sample.
/// @param fs Sampling frequency in samples/second.
/// @return Frequency in Hz.
double cycle_sample_2_hz(double f, double fs);

/// @brief Convert Hz to half cycles/sample.
/// @param f Frequency in cycles/second.
/// @param fs Sampling frequency in samples/second.
/// @return Half cycles/sample.
double hz_2_half_cycle_sample(double f, double fs);

/// @brief Convert half cycles/sample to Hz.
/// @param f Frequency in half cycles/sample.
/// @param fs Sampling frequency in samples/second.
/// @return Frequency in Hz.
double half_cycle_sample_2_hz(double f, double fs);

// Distance Errors
double l1_norm_mean(double *v1, double *v2, int size);
double l2_norm_mean(double *v1, double *v2, int size);
double mean_square(double *v1, double *v2, int size);

// Random Number Generation
float random_float(float min, float max);

// Signal Generation
double *generate_n_sines(double freqs[], size_t num, size_t samples_num, double fs);

#endif // UTILS_H
