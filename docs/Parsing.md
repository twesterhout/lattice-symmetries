# Parser

This document describes the syntax of the expressions that can be used to
construct `Expr` objects.

## Grammar

```
<coeff> ::= <signed-real> | "(" <signed-real> <plus-or-minus> <real> <imaginary-unit> ")"
<signed-real> ::= "-" <real> | <real>
<real> ::= <float> | <integer> "/" <integer> | <integer>
<imaginary-unit> ::= "im" | "j" | "ⅈ"
<plus-or-minus> ::= "+" | "-"
```
