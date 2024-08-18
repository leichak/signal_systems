#include <stdio.h>

#include "AnalogFilters.h"
#include "BillinearTransform.h"
#include "Filters.h"
#include "utils.h"

int main() {
  // // Steps generating test
  // const int step_num = 1000;
  // float start = 0.0;
  // float stop = 1.1;
  // float step = stop / (float)(step_num);
  // float steps[step_num] = {0.0};
  // fill_n_with_step(steps, start, step, step_num);

  // // Butterworth gain
  // float gains[step_num] = {0.0};

  // for (int i = 0; i < step_num; i++)
  // {
  //     gains[i] = butter_gain_low_pass_log10(3, 1.0, 0.5, steps[i]);
  //     printf("w rad/s: %f g: %f [dB]\n", steps[i], gains[i]);
  // }

  // Test calls
  test_coefficients_butter_sum_form_poles_coefficients();
  test_transform_wc();
  test_analog_to_digital();
}