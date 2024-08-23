#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <string.h>

#include "PbPlots/pbPlots.h"
#include "PbPlots/supportLib.h"

void HSLToRGB(double h, double s, double l, unsigned char *r, unsigned char *g, unsigned char *b);

void generate_colors(int m, RGBA *colors);

void *concat_path(char *path, char *prefix, char *filename);

void fill_n_with_step(double *xs, size_t size, float start, float stop);

int max_int(int a, int b);

ScatterPlotSeries *get_series_solid_thickness_color(double *xs, double *ys, size_t n, size_t thickness, RGBA *color);

int plot_x_y_single(double *xs, double *ys, size_t n, const char *prefix, char *filename);

int plot_x_y_overlay(double *xs[], double *ys[], size_t n_points, size_t overlayed_num, size_t thickness, char *labels[], const char *prefix, char *filename);

char *concat_strings(int count, ...);

wchar_t *char_to_wchar(const char *str);
#endif
