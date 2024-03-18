# Parser

This document describes the syntax of the expressions that can be used to
construct `Expr` objects.

## Grammar

### Complex numbers:

```
<coeff> ::= <signed-imaginary> | <signed-real> | "(" <signed-real> <plus-or-minus> <imaginary> ")"
<signed-imaginary> ::= "-" <imaginary> | <imaginary>
<signed-real> ::= "-" <real> | <real>
<imaginary> ::= <imaginary-unit> | <real> <imaginary-unit>
<real> ::= <float> | <integer> "/" <integer> | <integer>
<imaginary-unit> ::= "im" | "j" | "ⅈ"
<plus-or-minus> ::= "+" | "-"
```


