:bangbang: **Warning:** the following instructions are meant for developers.

## Reproducing

First, run the benchmarks:

```sh
python3 01_basis_construction.py
# or
# python3 02_operator_application.py
```

Next, plot the results:

```sh
latexrun 01_basis_construction.tex
# or
# python3 02_operator_application.tex
```

Note that these scripts expect data files to be in a specific folder. Modify
them for your own needs.

Finally, convert `.pdf` files to `.jpg` to include in the README:

```sh
pdftoppm -jpeg -r 600 01_basis_construction.pdf 01_basis_construction
mv 01_basis_construction-1.jpg 01_basis_construction.jpg
# or
# pdftoppm -jpeg -r 600 02_operator_application.pdf 02_operator_application
# mv 02_operator_application-1.jpg 02_operator_application.jpg
```
