#include "Timing.h"

void* c_timing_locale_data;

void* c_timing_get_locale_data(void) {
    return c_timing_locale_data;
}

void c_timing_set_locale_data(void* ptr) {
    c_timing_locale_data = ptr;
}
