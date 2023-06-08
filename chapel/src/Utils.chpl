module Utils {

use FFI;

proc initRuntime() {
  coforall loc in Locales do on loc {
    ls_hs_init();
  }
}

proc deinitRuntime() {
  coforall loc in Locales do on loc {
    ls_hs_exit();
  }
}


} // module Utils
