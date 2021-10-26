#include "lattice_symmetries_haskell.h"
#include <stdio.h>
#include <time.h>

int main(int argc, char **argv) {
  ls_hs_init();
  int x = ls_hs_foo(1);

  struct timespec t1;
  struct timespec t2;

  clock_gettime(CLOCK_REALTIME, &t1);
  for (int i = 0; i < 1000000; ++i) {
    x = ls_hs_foo(x);
  }
  clock_gettime(CLOCK_REALTIME, &t2);

  fprintf(stderr, "%i == %i\n", x, 1 + 123 + 1000000 * 123);

  unsigned long const delta =
      1.0e9 * (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec);
  fprintf(stderr, "Took: %lf ns\n", (double)delta / 1000000);

  ls_hs_exit();
  return 0;
}
