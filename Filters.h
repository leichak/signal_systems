// #include <stdarg.h>
// #include "Filters.c"
// #include "AnalogFilters.h"

// enum filter_type
// {
//     low_pass,
//     high_pass,
//     band_pass,
//     band_stop
// };

// enum analog_filter_design_type
// {
//     butterworth,
//     chebyshev_type1,
//     elliptic,
//     bessel
// };

// /// @brief This union represents various filters structs
// union filter
// {
//     struct
//     {
//         float cut_off
//     } low_pass;
//     struct
//     {
//         float cut_off
//     } high_pass;
//     struct
//     {
//         float cut_off_low;
//         float cut_off_high;
//     } band_pass;
//     struct
//     {
//         float cut_off_low;
//         float cut_off_high;
//     } band_stop;
// };

// /// @brief Function that takes design parameters and returns filter
// union filter get_discrete_filter_coefficients(enum filter_type t, float *cutoff[], float steepness, float ripple, float fs)
// {
//     switch (t)
//     {
//     case low_pass:
//         break;
//     case high_pass:
//         break;
//     case band_pass:
//         break;
//     case band_stop:
//         break;
//     }

//     union filter f;
//     return f;
// }
