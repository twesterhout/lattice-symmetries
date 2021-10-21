#include "lattice_symmetries_haskell.h"
#include <stdio.h>

int main() {
  ls_hs_init();
  printf("foo() = %i\n", ls_hs_foo());
  ls_hs_exit();
}
