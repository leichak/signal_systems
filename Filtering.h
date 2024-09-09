#ifndef FILTERING_H
#define FILTERING_H
/// @brief  A straightforward approach for IIR filter realization is direct form I,
/// where the difference equation is evaluated directly. This form is practical for small filters,
/// but may be inefficient and impractical (numerically unstable) for complex designs.[6] In general,
/// this form requires 2N delay elements (for both input and output signals) for a filter of order N.
/// Ref: https://en.wikipedia.org/wiki/Digital_filter
typedef struct
{
    float order
} DirectFormI;

typedef struct
{
    float order
} DirectFormII;

/// @brief Filtering using dirct form I
void process_df1(DirectFormI *ptr);

/// @brief  Filtering using direct form II
void process_df2(DirectFormII *ptr);

#endif // FILTERING_H
