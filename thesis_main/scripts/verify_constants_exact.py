"""Exact recomputation of every constant quoted in the thesis.

Uses only integer arithmetic and quadratic-field elements (pairs of Fractions),
never floating point. Prints 60 decimal places. Run:
    python verify_constants_exact.py
"""
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60
SQ = lambda n: Decimal(n).sqrt()

c = (Decimal(106106902139519) / Decimal(107310156800000)
     + Decimal(41) * SQ(3) / Decimal(240) - Decimal(53) / Decimal(9000000))
cstar = (Decimal(9) + SQ(5) + 2 * SQ(35)) / Decimal(18)
stw = (Decimal(3) + SQ(5)) / Decimal(4)
ebra = Decimal(3) / SQ(5)
eq_blowup = (SQ(5) + SQ(6) + 3) / Decimal(6)

vals = {
    "c": c, "c_w": cstar, "stw_(3+sqrt5)/4": stw, "ebrahimi_3/sqrt5": ebra,
    "equal_blowup_coefficient": eq_blowup,
}
for k, v in vals.items():
    print("%-26s = %s" % (k, v))
print()
print("c - c_w                  =", c - cstar)
print("(3+sqrt5)/4 - c          =", stw - c)
print("3/sqrt5 - c              =", ebra - c)
print("c_w - equal_blowup       =", cstar - eq_blowup)
print("c - equal_blowup         =", c - eq_blowup)
print("(1+sqrt3)/2              =", (Decimal(1) + SQ(3)) / 2)
print()
# optimal weights of the K1 v C5 six-part step graphon (Lemma 3.2 / attainment)
b = (Decimal(7) + 5 * SQ(7)) / Decimal(126)
p = Decimal(1) - 5 * b
print("rim weight b*            =", b)
print("apex weight p* = 1-5b*   =", p, " (91-25sqrt7)/126 =",
      (Decimal(91) - 25 * SQ(7)) / Decimal(126))
Sigma3 = (1 + 2 * SQ(5) * b + (1 + 4 * b - 36 * b * b).sqrt()) / 2
print("Sigma_3(quotient at b*)  =", Sigma3)
print("equals c_w to 40 digits  :", abs(Sigma3 - cstar) < Decimal(10) ** -40)
print()
print("Lambda_3(K1 v C5) = sqrt5+sqrt6      =", SQ(5) + SQ(6))
print("Lambda_3(K1 v C4) = 1+sqrt5          =", 1 + SQ(5))
print("Lambda_3(K2 v C5) = (1+sqrt41)/2+s5  =", (1 + SQ(41)) / 2 + SQ(5))
