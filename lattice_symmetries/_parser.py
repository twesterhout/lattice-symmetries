from lark import Lark, Transformer
import sympy
from sympy.physics.quantum import pauli
from sympy.physics.quantum.pauli import SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus

GRAMMAR = """
    ?start: sum_expr
    sum_expr: scaled_term*
    scaled_term: [SIGN] [coefficient] primitive_term+
    primitive_term: "(" sum_expr ")" | identity | spin | fermion
    identity: "I"
    spin: PREFIX SUPERSCRIPT subscript
    fermion: FERMION_OP subscript [SPIN_INDEX]
    coefficient: IMAGINARY | NUMBER | "(" NUMBER SIGN IMAGINARY ")"
    subscript: SUBSCRIPT_NUM | "_"? INT
    // Terminals
    PREFIX: "σ" | "S" | "\\\\sigma"
    SUPERSCRIPT: "ˣ" | "ʸ" | "ᶻ" | "⁺" | "⁻" | "^"? ("x"|"y"|"z"|"+"|"-")
    FERMION_OP: "c†" | "c^\\\\dagger" | "c" | "n"
    SPIN_INDEX: "↑" | "↓" | "_\\\\up" | "_\\\\down"
    SUBSCRIPT_NUM: /[₀₁₂₃₄₅₆₇₈₉]+/
    IMAGINARY: NUMBER ("j"|"ⅈ"|"im")
    SIGN: "+" | "-"
    
    %import common.INT
    %import common.NUMBER
    %import common.WS
    %ignore WS
"""

class ExpressionTransformer(Transformer):
    def sum_expr(self, items): return sympy.S.Zero if len(items) == 0 else sympy.Add(*items)
    def identity(self, items): return sympy.S.One

    def scaled_term(self, items):
        (sign, coeff, *terms) = items
        if coeff is None: coeff = sympy.S.One
        if sign == "-": coeff = -coeff
        return sympy.Mul(coeff, *terms)
    
    def primitive_term(self, items):
        if len(items) == 3:
            (prefix, expr, suffix) = items
            assert prefix == '(' and suffix == ')'
        else:
            (expr,) = items
        return expr

    def spin(self, items):
        (prefix, superscript, subscript) = items
        mapping = {"ˣ": SigmaX, "x": SigmaX, "ʸ": SigmaY, "y": SigmaY, "ᶻ": SigmaZ, "z": SigmaZ, "⁺": SigmaPlus, "+": SigmaPlus, "⁻": SigmaMinus, "-": SigmaMinus}
        assert isinstance(subscript, int)
        return sympy.Mul(sympy.Rational(1, 2) if prefix == "S" else sympy.S.One, mapping[superscript[-1]](subscript))

    def coefficient(self, items):
        if len(items) == 3:
            (real, sign, imag) = items
            real = sympy.Float(real.value)
            imag = sympy.Float(imag.value.rstrip('jⅈim')) * sympy.I
            sign = 1 if sign.value == "+" else -1
            return real + sign * imag
        elif items[0].type == "NUMBER":
            return sympy.Float(items[0].value)
        else:
            assert items[0].type == "IMAGINARY"
            value = items[0].value.rstrip('jⅈim')
            return sympy.Float(value) * sympy.I

    def subscript(self, items):
        if items[0].type == "SUBSCRIPT_NUM":
            mapping = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9"}
            return int("".join(map(mapping.__getitem__, items[0].value)), base=10)
        else:
            assert len(items) == 1 and items[0].type == "INT"
            return int(items[0].value, base=10)

_PARSER = Lark(GRAMMAR, parser='lalr', transformer=ExpressionTransformer())

def parse_expr(text: str): return _PARSER.parse(text)
