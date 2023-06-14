#include "lattice_symmetries_types.h"
#include <HsFFI.h>
#include <Rts.h>
// #include <stdio.h>

void ls_hs_init(void) {
  int argc = 1;
  char *argv[] = {"lattice_symmetries", NULL};
  // "+RTS", "-N1", "--install-signal-handlers=no", "-RTS", NULL};
  char **pargv = argv;

  // For some reason, options from argv are not processed properly, so we
  // manually set all RTS options using rts_opts field of RtsConfig
  RtsConfig conf = defaultRtsConfig;
  conf.rts_opts_enabled = RtsOptsAll;
  conf.rts_opts = "-N1 --install-signal-handlers=no";
  hs_init_ghc(&argc, &pargv, conf);
}

void ls_hs_exit(void) {
  hs_exit();
}
