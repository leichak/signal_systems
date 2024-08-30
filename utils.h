#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "PbPlots/pbPlots.h"
#include "PbPlots/supportLib.h"

#define TWICE_PI 2.0 * M_PI;
#define TWICE_PI_INVERSE 1.0 / 2.0 * M_PI;

// Strings utils
void *concat_path(char *path, char *prefix, char *filename);

char *concat_strings(int count, ...);

wchar_t *char_to_wchar(const char *str);

// Drawing utils
void HSLToRGB(double h, double s, double l, unsigned char *r, unsigned char *g, unsigned char *b);

ScatterPlotSeries *get_series_solid_thickness_color(double *xs, double *ys, size_t n, size_t thickness, RGBA *color);

int plot_x_y_single(double *xs, double *ys, size_t n, const char *prefix, char *filename);

int plot_x_y_overlay(double *xs[], double *ys[], size_t n_points, size_t overlayed_num, size_t thickness, char *labels[], const char *prefix, char *filename);

void generate_colors(int m, RGBA *colors);

// Pre computed utils

void fill_n_with_step(double *xs, size_t size, float start, float stop);

// Numerical utils

int max_int(int a, int b);

// Units utils

/// @brief Function converting rad/s to hz. 1 Hz equals 2PI rad/s
/// @param omega omega in PI radians / s
/// @return Hz - cycles / s
double rad_s_2_hz(double omega);

/// @brief Function converting hz to rad/s. 1 Hz equals 2PI rad/s
/// @param f Hz - cycles / s
/// @return omega in PI radians / s
double rad_s_2_hz_inv(double f);

/// @brief Function converting rad/s to rad/sample
/// @param omega omega in PI radians / s
/// @param fs sampling frequency in sample / s
/// @return radians / sample
double rad_s_2_rad_sample(double omega, double fs);

/// @brief Function converting rad/sample to rad/s
/// @param omega radians / sample
/// @param fs sampling frequency in sample / s
/// @return PI radians / s
double rad_s_2_rad_sample_inv(double omega, double fs);

/// @brief Function converting Hz to cycle/sample
/// @param f frequency cycle / second
/// @param fs sampling frequency sample / second
/// @return cycle / sample
double hz_2_cycle_sample(double f, double fs);

/// @brief Function converting cycle/sample to Hz
/// @param f frequency cycle / sample
/// @param fs sampling frequency sample / second
/// @return Hz
double hz_2_cycle_sample_inv(double f, double fs);

/// @brief Function converting Hz to half cycle/sample
/// @param f frequency cycle / second
/// @param fs sampling frequency sample / second
/// @return half cycle / sample
double hz_2_half_cycle_sample(double f, double fs);

/// @brief Function converting half cycle/sample to Hz
/// @param f frequency half cycle / sample
/// @param fs sampling frequency sample / second
/// @return Hz
double hz_2_half_cycle_sample_inv(double f, double fs);

#endif
