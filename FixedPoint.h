
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP

#include <math.h>
/*
Source:

Increased complexity of floating point:
    1. Power consumption of target processor if this supports floating point.
    2. If the program needs to be as fast the energy will increase.

We convert floating point the nto fixed point data representation and call it
fixed point refinement.

Unsigned fixed point representation:
    Unsigned integer with bits bn-1 bn-2 ... b1 b0
    Value of unsigned integer:
        Vuint<n> = bn-1 * 2^n + bn-2 * 2^n-2 + .. + b1 * 2 + b0 * 1

In fixed point representation binary point shifts from right most position to a different position
k.
Fixed point representation is defined by the number of bits N as well as the position of binary
point k. For an unsigned number we adopt the notation ufix<N,k>.

Integer defined above would be ufix<N,0>.

Value of ufix<N,k>:
    Vufix<N,k> = bN-k-1 * 2^N-k-1 + bN-k-2 * 2^N-k-2 + ... + b0 + b-1*2^-1 + .. + bk2^-k

An ufix<N,k> fits in the same number of bits as N-bit unsigned integer, only difference in our
interpretation of what these bits mean. Following we have an example of ufix<6,4>. The given
bit pattern can be evaluated to the formulas above to compute the integer value as well as the
fixed point value.

2^1 2^0 2^-1 2^-2 2^-3 2^-4
b1  b0  b-1  b-2  b-3  b-4  | unsigned  |  ufix<6,4>
                            |           |
0   0   0    0    0    0    |    0      |     0
0   0   0    0    0    1    |    1      |     0.0625
0   0   1    0    0    1    |    9      |     0.5625
0   1   0    0    0    0    |   16      |     1.
1   1   1    1    1    1    |   63      |     3.9375

There is a relationship between the value of unsigned integer of N bits and the value of unsigned
fixed point number ufix<N,k>:

    Vuint<N> = Vufix<N,k>*2^k
    Vufix<N,k> = Vuint<N>*2^-k

The floating point number is converted to a fixed-point number by proper scaling and conversion to
an integer data type.

Signed Fixed-point representation

Signed values in integer arithmetic are captured in one of three different ways, the most common
method is two's complement and we will focus on this method. 2 others methods are one's complement
and sign-magnitude representation.

In a two complements representation the value of signed integer is defined as follows.

    Vint<N> = bN-1*[-(2^N-1)] + bN-2*2^N-2 + .. + b1*2 + b0

The weight of most significant bit is negative, and fixed point  has it similarly:

    Vfix<N,k> bN-k-1[-(2^N-k-1)] + bN-k-2*2^N-2 + ... + b0 + b-1 * 2^-1 + .. + bk 2^-k

2^1 2^0 2^-1 2^-2 2^-3 2^-4
b1  b0  b-1  b-2  b-3  b-4  |   signed  |  ufix<6,4>
                            |           |
0   0   0    0    0    0    |    0      |     0
0   0   0    0    0    1    |    1      |     0.0625
1   0   1    0    0    0    |   -24     |     -1.5
1   1   0    0    0    0    |   -16     |     -1.0
1   1   1    1    1    1    |   -1      |     -0.0625

Relationships between integer value and fixed point representation.

    Vint<N> = Vfix<N,k>*2^k
    Vfix<N,k> = Vint<N>*2^-k

1. First scale it by 2^k so that the least significant fractional bit gets unit weight. This
conversion can lead to precision loss (at the LSB side) and overflow (at the MSB side).

Example:

#include <stdio.h>

    float taps[5] = {0.1, 1.2, -0.3, 3.4, -1.5};
    int taps_8_6[5]; // a fix <8,6> type
    unsigned i;
    for (i=0; i<5; i++)
        taps_8_6[i] (int) (taps[i] * (1 << 6)); // multiply by 2^k

    for (i=0; i<5; i++)
        printf("%+4.2f %10d\n", taps[i], taps_8_6[i]);

    Output:
    +0.10   6
    +1.20   76
    -0.30   -19
    +3.40   217
    -1.50   -96

The value of fix<8,6> is between -2 10000000 and 1.98... for 011111111, hence the floating value
+3.40 is overflow.

Some of the floating point numbers like 0.1 cannot be expressed with full precision, which means
for example that value 6 maps to 00000110 or 0.09...

fix<8,6> - means 8 bits total with 6 fractional bits
Scale it by 2^k, which means that fractional part will be represented as integer and lsb
will have unit weight - each bit represents unit of two.

Sign extension is hte process of increasing the number of bits used to represent a number
while preserving its sign and value. When a fixed point number represented in two's complement
is converted to a larger integer format (8 to 32 bit) the sign bit is replicated into the new bits
to maintain the correct negative value. Otherwise it could be misinterpreted.

So basically we convert out floating point to fixed using the least N bits with k fractional bits
And if we use less bits than our N then for signed we have sign extension.

After computations in fixed point precision we want to convert back to floating point.
So the least significant bit gets it actual weight.


If location of MSB of the fixed point representation is not located at the MSB integer
representation we need to replicated MSB of fixed point to integer representation to get sign
correctly like:

11111111 11111111 11111111 11000011, last 8 bits is fixed point and 1 bit of those 8 is
replicated across all bits

Example code:

void fn fixed_2_float() {
    int taps_8_6[5] = {6, 75, -19, 127, -96};
    unsigned i;

    for (i=0; i<5; i++)
        taps[i] = (taps_8_6[i] * (1.0 / (1 << 6))); // multiply by 2^-kk

    for (i=0; i<5; i++)
        printf("%+4.2f %10d\n", taps[i], taps_8_6[i]);
}

Fixed point arithmetic - unsigned numbers

Unsigned addition:

When we add two ufix<N,k> numbers then the result is an ufix<N+1,k> number, the extra bit on
MSB side is to capture carry bit in case one is generated.
The subtraction of 2 numbers uses the same rule since subtraction can be defined as addition of a
ufix<N,k> with the two complement version of the other ufix<N,k>

We need to align ufix<N,k1> and ufix<N,k2> k1 > k2 before those can be added, this will increase wordlength
of the sum with k1 - k2 bits.

If we want to capture the sum as N+1 bit number then there are two possible alignments,
1. For ufix<N+1,k1> sum we will have to increase the number of ufix<N,k2> fractional bits with
a left shift before addition
2. For ufix<N+1,k2> sum we will have to decrease the the number of ufix<N,k1> fractional bits with a
right shift before addition

Since total sum has only N + 1 bits there is potential precision loss:
    - left shifting V2 may case an overflow when MSB side bits are lost
    - while right shifting V1 may cause precision loss when LSB-side bits are lost
    - The bottom line is that when combining numbers with a different number of fraciotnal bits
    an alignment must be done, which can cause overflow or precision loss


Multiplication - when we multiply two ufix<N,k> numbers then the result is a ufix<2N,2k>
If the result has to be captured in an ufix<N,k> then result must be right shifted by k bits.
Precision loss may occur both at MSB side as well at the LSB side

When computing a digital filter using fixed point arithmetic, the precision loss is more pronounced
than with floating point arithmetic. Errors add to quantization noise. It is uniform distributed
random variable so it appears like a broadband noise in useful signal.
*/

int to_fix(float x, int k)
{
    return (int)(x * (1 << k)); // equivalent 2^k
}

float to_float(float x, int k)
{
    return (x * 1.0 / (1 << k)); // equivalent 2^-k
}

int fixed_N_k_mul(int x, int y, int k) // result is 2N
{
    return (x * y) >> k;
}

int fixed_N1_k_add(int x, int y) // N+1 k result
{
    return (x + y);
}

int fixed_N1_k2_add(int x_k1, int y_k2, int k1, int k2) // loss of precision
{                                                       // k1 > k2
    return ((x_k1 >> abs(k1 - k2)) + y_k2);
}

int fixed_N1_k1_add(int x_k1, int y_k2, int k1, int k2) // overflow
{                                                       // k1 > k2
    return (x_k1 + (y_k2 << (abs(k1 - k2))));
}

#endif
