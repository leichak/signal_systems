#include <stdio.h>

#include "Tests.h"

int main()
{
    // Tests
    test_generate_filters();
    test_various_orders_filters();
    test_analog_to_digital();
    test_make_causal();
}
