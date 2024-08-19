#include <stdio.h>

#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "utils.h"

int main()
{
    // Tests
    test_generate_filters();
    test_various_orders_filters();
    test_analog_to_digital();
}
