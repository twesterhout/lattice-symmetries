This folder contains files which were used to measure performance of `SPINPACK`.

* Version: spinpack-2.58a;
* [`Make_cfg.inc`](./Make_cfg.inc) specifies the compilation options (important
  options are `-march=native`, `-DCONFIG_PTHREAD=64`, and `-DVecType=1`);
* `./spinpack/exe/m_square 6 0 0 6 1 > daten.def` command was used to generate
  [`daten.def`](./daten.def) input file.
* `./spinpack/exe/spin -t64 --maxmem=0` command was used to perform the
  diagonalization.
* Timings are taken from `tmp/tri.txt` using the following command:
  ~~~~~~sh
  cat tmp/tri.txt | sed -E 's/\s+/,/g' | cut -d',' -f 5 | tr -d m > tri.txt.processed
  ~~~~~~
* We then fit a line through it to obtain `0.9244` minutes per iterations, i.e.
  `55.5` seconds.
