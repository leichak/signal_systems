
#include "DigitalFilters.h"

void free_digital_filter(DigitalFilter *p) {
  if (p != NULL) {
    free(p->a_k);
    free(p->b_k);
  }
}