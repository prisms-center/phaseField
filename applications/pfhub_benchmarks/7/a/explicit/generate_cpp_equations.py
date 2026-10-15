#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 12:47:25 2026

@author: staublin
"""

import sympy
from sympy import cos, sech, sin, sqrt, tanh

x, y = sympy.symbols("x y", real=True)
t = sympy.symbols("t", positive=True, real=True)
kappa = sympy.symbols("kappa", positive=True, real=True)
A_1, A_2 = sympy.symbols("A_1 A_2", real=True)
B_1, B_2 = sympy.symbols("B_1 B_2", real=True)
C_2 = sympy.symbols("C_2", real=True)


# ------------- These Lines Copied From PFHub Benchmark 7 -------------------
S = (
    -(
        2
        * sqrt(kappa)
        * (A_1 * B_1 * t * cos(B_1 * x) + A_2 * B_2 * cos(B_2 * x + C_2 * t)) ** 2
        * tanh(
            sqrt(2)
            * (-A_1 * t * sin(B_1 * x) - A_2 * sin(B_2 * x + C_2 * t) + y - 1 / 4)
            / (2 * sqrt(kappa))
        )
        + sqrt(2)
        * kappa
        * (-A_1 * B_1**2 * t * sin(B_1 * x) - A_2 * B_2**2 * sin(B_2 * x + C_2 * t))
        - sqrt(2) * (A_1 * sin(B_1 * x) + A_2 * C_2 * cos(B_2 * x + C_2 * t))
    )
    * sech(
        sqrt(2)
        * (-A_1 * t * sin(B_1 * x) - A_2 * sin(B_2 * x + C_2 * t) + y - 1 / 4)
        / (2 * sqrt(kappa))
    )
    ** 2
    / (4 * sqrt(kappa))
)
alpha = A_1 * t * sin(B_1 * x) + A_2 * sin(B_2 * x + C_2 * t) + 1 / 4
eta = (
    tanh(
        sqrt(2)
        * (4 * A_1 * t * sin(B_1 * x) + 4 * A_2 * sin(B_2 * x + C_2 * t) - 4 * y + 1)
        / (8 * sqrt(kappa))
    )
    / 2
    + 1 / 2
)
eta0 = (
    tanh(sqrt(2) * (4 * A_2 * sin(B_2 * x) - 4 * y + 1) / (8 * sqrt(kappa))) / 2 + 1 / 2
)
# ------------- These Lines Copied From PFHub Benchmark 7 -------------------

print("\neta0 = ")
print(sympy.cxxcode(eta0))

print("\neta = ")
print(sympy.cxxcode(eta))

print("\nalpha = ")
print(sympy.cxxcode(alpha))

# Handle S differently because sympy likes to interpret sech(...) as complex
# exponentials instead of just doing 1 / cosh(...)
print("\nS = ")
s_str = str(S)
s_str = s_str.replace("sqrt(", "std::sqrt(")
s_str = s_str.replace("sin(", "std::sin(")
s_str = s_str.replace("cos(", "std::cos(")
s_str = s_str.replace("tanh(", "std::tanh(")
print(s_str)
