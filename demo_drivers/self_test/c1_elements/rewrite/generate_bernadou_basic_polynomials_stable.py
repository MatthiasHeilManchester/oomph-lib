#!/usr/bin/env python3
"""Generate stable, loop-free Bernadou basic-polynomial evaluators.

Usage:
  python3 generate_bernadou_basic_polynomials.py matrix_3.cc 3 7 > basic_3.h
  python3 generate_bernadou_basic_polynomials.py matrix_5.cc 5 9 > basic_5.h

The input file must contain exactly one monomial_to_basic_matrix body, using
assignments A(i,j)=...; or b2l(i,j)=...;. Coefficients may be rational numbers
and rational multiples of sqrt(2).

For every value/derivative component this program:
  1. forms the polynomial exactly in Q(sqrt(2))[s,t];
  2. extracts all common powers of s, t, and u=1-s-t;
  3. evaluates the residual by bivariate Horner nesting;
  4. chooses the cheaper of s-outer and t-outer Horner forms; and
  5. expands the result exactly and verifies it against the original.
"""

import argparse
import ast
import re
import sys
from fractions import Fraction

ASSIGNMENT = re.compile(
    r"\b(?:A|b2l)\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*=\s*(.*?);",
    re.DOTALL)


class Q2:
    """Exact a+b*sqrt(2), with rational a and b."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, other):
        other = as_q2(other)
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-as_q2(other))

    def __rsub__(self, other):
        return as_q2(other) - self

    def __mul__(self, other):
        other = as_q2(other)
        return Q2(self.a * other.a + 2 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_q2(other)
        denominator = other.a * other.a - 2 * other.b * other.b
        if denominator == 0:
            raise ZeroDivisionError
        return Q2((self.a * other.a - 2 * self.b * other.b) / denominator,
                  (self.b * other.a - self.a * other.b) / denominator)

    def __rtruediv__(self, other):
        return as_q2(other) / self

    def __eq__(self, other):
        other = as_q2(other)
        return self.a == other.a and self.b == other.b

    def is_zero(self):
        return self.a == 0 and self.b == 0


def as_q2(value):
    return value if isinstance(value, Q2) else Q2(value)


def parse_coefficient(text):
    replacements = (
        ("std::sqrt(2.0)", "SQRT_TWO"),
        ("std::sqrt(2.)", "SQRT_TWO"),
        ("std::sqrt(2)", "SQRT_TWO"),
        ("sqrt(2.0)", "SQRT_TWO"),
        ("sqrt(2.)", "SQRT_TWO"),
        ("sqrt(2)", "SQRT_TWO"),
    )
    for old, new in replacements:
        text = text.replace(old, new)
    node = ast.parse(text.strip(), mode="eval").body

    def evaluate(item):
        if isinstance(item, ast.Constant) and isinstance(item.value, (int, float)):
            return Q2(Fraction(str(item.value)))
        if isinstance(item, ast.Name) and item.id == "SQRT_TWO":
            return Q2(0, 1)
        if isinstance(item, ast.UnaryOp) and isinstance(item.op, ast.USub):
            return -evaluate(item.operand)
        if isinstance(item, ast.UnaryOp) and isinstance(item.op, ast.UAdd):
            return evaluate(item.operand)
        if isinstance(item, ast.BinOp):
            left = evaluate(item.left)
            right = evaluate(item.right)
            if isinstance(item.op, ast.Add):
                return left + right
            if isinstance(item.op, ast.Sub):
                return left - right
            if isinstance(item.op, ast.Mult):
                return left * right
            if isinstance(item.op, ast.Div):
                return left / right
        raise ValueError("Unsupported coefficient expression: " + ast.dump(item))

    return evaluate(node)


# A polynomial is a dictionary {(power_of_s, power_of_t): exact_coefficient}.
def clean(poly):
    return {key: value for key, value in poly.items() if not value.is_zero()}


def add_poly(left, right):
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, Q2()) + value
    return clean(result)


def multiply_poly(left, right):
    result = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            key = (a + c, b + d)
            result[key] = result.get(key, Q2()) + x * y
    return clean(result)


def shift_poly(poly, ds, dt):
    return clean({(a + ds, b + dt): value for (a, b), value in poly.items()})


def differentiate(poly, ds, dt):
    result = {}
    for (a, b), value in poly.items():
        if a < ds or b < dt:
            continue
        multiplier = 1
        for k in range(ds):
            multiplier *= a - k
        for k in range(dt):
            multiplier *= b - k
        result[(a - ds, b - dt)] = value * multiplier
    return clean(result)


def monomial_powers(order):
    return [(degree - t_power, t_power)
            for degree in range(order, -1, -1)
            for t_power in range(degree + 1)]


def extract_st_factors(poly):
    if not poly:
        return 0, 0, {}
    s_power = min(a for a, unused in poly)
    t_power = min(b for unused, b in poly)
    residual = {(a - s_power, b - t_power): value
                for (a, b), value in poly.items()}
    return s_power, t_power, clean(residual)


U_POLY = {(0, 0): Q2(1), (1, 0): Q2(-1), (0, 1): Q2(-1)}


def divide_once_by_u(poly):
    """Return exact quotient by 1-s-t, or None if not divisible."""
    if not poly:
        return {}

    # Regard p as a polynomial in t. If q=sum q_b(s)t^b, then
    # p=(1-s-t)q and, descending in b,
    # q_{b-1}(s)=(1-s)q_b(s)-p_b(s).
    max_t = max(b for unused, b in poly)
    p_by_t = {}
    for (a, b), value in poly.items():
        p_by_t.setdefault(b, {})[a] = value

    def add_s_poly(x, y):
        result = dict(x)
        for power, value in y.items():
            result[power] = result.get(power, Q2()) + value
        return {k: v for k, v in result.items() if not v.is_zero()}

    def minus_s_times(x):
        # (1-s)*x
        result = dict(x)
        for power, value in x.items():
            result[power + 1] = result.get(power + 1, Q2()) - value
        return {k: v for k, v in result.items() if not v.is_zero()}

    q_by_t = {}
    q_b = {}
    for b in range(max_t, 0, -1):
        p_b = p_by_t.get(b, {})
        q_previous = add_s_poly(minus_s_times(q_b),
                                {k: -v for k, v in p_b.items()})
        q_by_t[b - 1] = q_previous
        q_b = q_previous

    # The constant coefficient must satisfy p_0=(1-s)q_0.
    remainder = add_s_poly(p_by_t.get(0, {}),
                           {k: -v for k, v in minus_s_times(q_by_t.get(0, {})).items()})
    if remainder:
        return None

    quotient = {}
    for b, s_polynomial in q_by_t.items():
        for a, value in s_polynomial.items():
            quotient[(a, b)] = value
    quotient = clean(quotient)

    # Defensive exact reconstruction check for the division itself.
    if multiply_poly(quotient, U_POLY) != clean(poly):
        raise RuntimeError("Internal error in exact division by 1-s-t")
    return quotient


def extract_u_factors(poly):
    power = 0
    residual = clean(poly)
    while residual:
        quotient = divide_once_by_u(residual)
        if quotient is None:
            break
        residual = quotient
        power += 1
    return power, residual


def number(text):
    if text.denominator == 1:
        return f"{text.numerator}.0"
    return f"({text.numerator}.0 / {text.denominator}.0)"


def coefficient_cpp(value):
    if value.b == 0:
        return number(value.a)
    if value.a == 0:
        if value.b == 1:
            return "sqrt_two"
        if value.b == -1:
            return "-sqrt_two"
        return f"{number(value.b)} * sqrt_two"
    sign = "+" if value.b > 0 else "-"
    return f"({number(value.a)} {sign} {number(abs(value.b))} * sqrt_two)"


class Expression:
    __slots__ = ("text", "multiplies", "adds")

    def __init__(self, text, multiplies=0, adds=0):
        self.text = text
        self.multiplies = multiplies
        self.adds = adds

    @property
    def cost(self):
        return self.multiplies + self.adds


def constant_expression(value):
    return Expression(coefficient_cpp(value))


def horner_univariate(coefficients, variable):
    highest = max(coefficients, default=0)
    expression = constant_expression(coefficients.get(highest, Q2()))
    for power in range(highest - 1, -1, -1):
        coefficient = constant_expression(coefficients.get(power, Q2()))
        expression = Expression(
            f"({coefficient.text} + {variable} * {expression.text})",
            coefficient.multiplies + expression.multiplies + 1,
            coefficient.adds + expression.adds + 1)
    return expression


def horner_bivariate(poly, outer):
    """Hornerise in inner variable, then in outer variable."""
    if not poly:
        return Expression("0.0")
    if outer == "t":
        outer_index, inner_index, inner_name = 1, 0, "s"
    else:
        outer_index, inner_index, inner_name = 0, 1, "t"

    grouped = {}
    for powers, value in poly.items():
        outer_power = powers[outer_index]
        inner_power = powers[inner_index]
        grouped.setdefault(outer_power, {})[inner_power] = value

    highest = max(grouped)
    expression = horner_univariate(grouped.get(highest, {}), inner_name)
    for power in range(highest - 1, -1, -1):
        coefficient = horner_univariate(grouped.get(power, {}), inner_name)
        expression = Expression(
            f"({coefficient.text} + {outer} * {expression.text})",
            coefficient.multiplies + expression.multiplies + 1,
            coefficient.adds + expression.adds + 1)
    return expression


def power_name(variable, power):
    if power == 0:
        return None
    if power == 1:
        return variable
    return f"{variable}{power}"


def factorised_expression(poly):
    if not poly:
        return Expression("0.0"), (0, 0, 0), {}

    s_power, t_power, residual = extract_st_factors(poly)
    u_power, residual = extract_u_factors(residual)

    # Removing u factors can reveal no new common s/t factor: u has a nonzero
    # constant term. The order above is therefore sufficient and canonical.
    by_t = horner_bivariate(residual, "t")
    by_s = horner_bivariate(residual, "s")
    expression = min((by_t, by_s), key=lambda item: (item.cost, len(item.text)))

    factors = [power_name("s", s_power), power_name("t", t_power),
               power_name("u", u_power)]
    factors = [item for item in factors if item is not None]
    text = expression.text
    multiplies = expression.multiplies
    for factor in reversed(factors):
        text = f"{factor} * ({text})"
        multiplies += 1

    # Exact reconstruction, independent of the emitted Horner syntax.
    reconstructed = residual
    reconstructed = shift_poly(reconstructed, s_power, t_power)
    for unused in range(u_power):
        reconstructed = multiply_poly(reconstructed, U_POLY)
    if reconstructed != clean(poly):
        raise RuntimeError("Exact factorisation reconstruction failed")

    return Expression(text, multiplies, expression.adds), (s_power, t_power, u_power), residual


def wrap_assignment(lhs, expression, width=92):
    prefix = f"  {lhs} = "
    if len(prefix) + len(expression.text) <= width:
        return [prefix + expression.text + ";"]
    # The nested expression remains mechanically generated; wrapping at
    # operators is deliberately left to clang-format.
    return [f"  {lhs} =", f"    {expression.text};"]


def generate_function(boundary_order, name, output_type, output_name,
                      derivatives, rows):
    expressions = []
    maximum_s = maximum_t = maximum_u = 0
    uses_sqrt_two = False
    total_cost = 0

    for row in rows:
        for derivative in derivatives:
            poly = differentiate(row, *derivative)
            expression, factors, residual = factorised_expression(poly)
            expressions.append(expression)
            maximum_s = max(maximum_s, factors[0])
            maximum_t = max(maximum_t, factors[1])
            maximum_u = max(maximum_u, factors[2])
            uses_sqrt_two |= any(value.b != 0 for value in poly.values())
            total_cost += expression.cost

    lines = [
        "template<>",
        f"void BernadouElementBasis<{boundary_order}>::{name}(",
        f"  const Vector<double>& s_basic, {output_type}& {output_name}) const",
        "{",
        "  const double s = s_basic[0];",
        "  const double t = s_basic[1];",
    ]
    if maximum_u:
        lines.append("  const double u = 1.0 - s - t;")

    for variable, maximum in (("s", maximum_s), ("t", maximum_t), ("u", maximum_u)):
        for power in range(2, maximum + 1):
            previous = variable if power == 2 else f"{variable}{power - 1}"
            lines.append(f"  const double {variable}{power} = {previous} * {variable};")
    if uses_sqrt_two:
        lines.append("  static const double sqrt_two = std::sqrt(2.0);")
    lines.append("")

    index = 0
    for row_number in range(len(rows)):
        for column in range(len(derivatives)):
            lhs = (f"{output_name}[{row_number}]" if len(derivatives) == 1
                   else f"{output_name}({row_number}, {column})")
            lines.extend(wrap_assignment(lhs, expressions[index]))
            index += 1
    lines.extend(["}", ""])
    return lines, total_cost


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("boundary_order", type=int)
    parser.add_argument("basis_order", type=int)
    args = parser.parse_args()

    text = open(args.source, encoding="utf-8").read()
    matches = ASSIGNMENT.findall(text)
    if not matches:
        raise SystemExit("No A(i,j)=... or b2l(i,j)=... assignments were found")

    dimensions = (args.basis_order + 1) * (args.basis_order + 2) // 2
    entries = {}
    for row_text, column_text, coefficient_text in matches:
        key = (int(row_text), int(column_text))
        if key in entries:
            raise SystemExit(f"Duplicate matrix assignment for entry {key}")
        if key[0] >= dimensions or key[1] >= dimensions:
            raise SystemExit(f"Matrix entry {key} lies outside {dimensions}x{dimensions}")
        entries[key] = parse_coefficient(coefficient_text)

    powers = monomial_powers(args.basis_order)
    rows = []
    for row in range(dimensions):
        polynomial = {}
        for column, power in enumerate(powers):
            value = entries.get((row, column), Q2())
            if not value.is_zero():
                polynomial[power] = value
        rows.append(clean(polynomial))

    output = [
        "// Automatically generated from monomial_to_basic_matrix().",
        "// Coefficients and all symbolic transformations were exact in Q(sqrt(2)).",
        "// Common s, t and (1-s-t) factors were extracted before Horner evaluation.",
        "",
    ]
    specifications = (
        ("full_basic_polynomials", "Shape", "phi", ((0, 0),)),
        ("dfull_basic_polynomials", "DShape", "dphi", ((1, 0), (0, 1))),
        ("d2full_basic_polynomials", "DShape", "d2phi",
         ((2, 0), (1, 1), (0, 2))),
    )

    for name, output_type, output_name, derivatives in specifications:
        generated, cost = generate_function(
            args.boundary_order, name, output_type, output_name,
            derivatives, rows)
        output.extend(generated)
        print(f"Verified {name}: exact reconstruction passed; "
              f"estimated {cost} additions/multiplications.", file=sys.stderr)

    print("\n".join(output))


if __name__ == "__main__":
    main()
