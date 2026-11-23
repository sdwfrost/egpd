#!/usr/bin/env python3
"""
Third/fourth derivatives (d34) of the continuous EGPD negative log-likelihood.

WHY THIS EXISTS
---------------
`egpd1d34`..`egpd4d34` in src/egpd.cpp were generated long ago by another tool;
egpd1d34 and egpd3d34 feed `.d1H0()` (R/logdetH.R), i.e. the analytic REML
gradient used by the default outer optimiser. Their entry (1,2,2) =
d^3 NLL / (d lpsi d xi d xi) was wrong (round-3 review finding N1/E-02): e.g.
carrier 1 at lpsi = log 2, xi = 0.3, lkappa = log 1.5, y = 3.3 gave -1.30501
instead of -1.08032, so the REML gradient was not the gradient of the REML
criterion. SymPy shows the error is in the GPD-density part, common to every
carrier: with z = y/psi,
    d^3 [lpsi + (1 + xi) A] / (d lpsi d xi^2) = 2 z^2 (1 - z) / (1 + xi z)^3,
whereas the generated code had 4 z (1 - z) / (1 + xi z)^3. That term was
corrected in place in all four functions. egpd4d34 also had a wrong fourth
derivative (2,2,3,3); it is replaced by the helper emitted by
`entry egpd4 2,2,3,3`. Every other entry agrees with the symbolic values
(mode `check`) to ~1e-12 and is unchanged.

PARAMETERS AND OBJECTIVES (y > 0, 1 + xi*y/psi > 0)
---------------------------------------------------
  A(y)  = -log(1 - H(y)) = log1p(xi*y/psi)/xi          (GPD tail exponent)
  -log h(y) = lpsi + (1 + xi)*A                         (GPD density part)
  carrier 1 (lpsi, xi, lkappa):          G(v) = v^kappa
      NLL = lpsi + (1+xi)A - lkappa - (kappa-1) log(1 - e^{-A})
  carrier 3 (lpsi, xi, ldelta):          g(v) = (1 + 1/delta)(1 - (1-v)^delta)
      NLL = lpsi + (1+xi)A - log1p(1/delta) - log(1 - e^{-delta A})
  carrier 4 (lpsi, xi, ldelta, lkappa):  G(v) = G3(v)^(kappa/2)
  carrier 2, direct parametrisation used only by egpd2d34
            (lpsi, xi, lkappa1, lkappa2, logitp): G(v) = p v^k1 + (1-p) v^k2

COLUMN ORDER (must match .indices() in R/logdetH.R)
--------------------------------------------------
Graded lexicographic order of the derivative multi-index, third derivatives
first; for 3 parameters:
  d3:  (3,0,0) (2,1,0) (2,0,1) (1,2,0) (1,1,1) (1,0,2) (0,3,0) (0,2,1) (0,1,2) (0,0,3)
i.e. (111,112,113,122,123,133,222,223,233,333); then the 15 fourth derivatives.

USAGE
-----
    python derive_egpd_d34.py entry egpd1 1,2,2     # C++ helper for one entry
    python derive_egpd_d34.py entry egpd4 2,2,3,3
    python derive_egpd_d34.py gpd                   # the GPD-density (1,2,2) term
    python derive_egpd_d34.py check                 # 40-digit reference values

The `check` output is the source of the hard-coded reference constants in
tests/testthat/test-audit-native-r3.R.
"""

import sys
import itertools
import sympy as sp
from sympy import symbols, exp, log
from sympy.codegen.cfunctions import log1p, expm1
try:                                     # sympy >= 1.8
    from sympy.printing.c import ccode
except ImportError:                      # older sympy
    from sympy.printing.ccode import ccode

y = symbols('y', positive=True)
lpsi, xi = symbols('lpsi xi', real=True)
t2, t3, t4 = symbols('t2 t3 t4', real=True)


def tail_exponent():
    return log1p(xi * y * exp(-lpsi)) / xi


def egpd1_nll():
    A = tail_exponent()
    k = exp(t2)
    return lpsi + (1 + xi) * A - t2 - (k - 1) * log(-expm1(-A))


def egpd3_nll():
    A = tail_exponent()
    d = exp(t2)
    return lpsi + (1 + xi) * A - log1p(exp(-t2)) - log(-expm1(-d * A))


def egpd4_nll():
    A = tail_exponent()
    d, k = exp(t2), exp(t3)
    G3 = 1 - exp(-A) * (1 + (1 - exp(-d * A)) / d)
    return (lpsi + (1 + xi) * A - log1p(exp(-t2)) - log(-expm1(-d * A))
            - (t3 - log(2)) - (k / 2 - 1) * log(G3))


def egpd2_direct_nll():
    A = tail_exponent()
    k1, k2 = exp(t2), exp(t3)
    p = 1 / (1 + exp(-t4))
    H = -expm1(-A)
    g = p * k1 * H ** (k1 - 1) + (1 - p) * k2 * H ** (k2 - 1)
    return lpsi + (1 + xi) * A - log(g)


MODELS = {
    'egpd1': (egpd1_nll, (lpsi, xi, t2)),
    'egpd3': (egpd3_nll, (lpsi, xi, t2)),
    'egpd4': (egpd4_nll, (lpsi, xi, t2, t3)),
    'egpd2direct': (egpd2_direct_nll, (lpsi, xi, t2, t3, t4)),
}

# C++ argument names of the carrier parameters (t2, t3, ...).
PAR_NAMES = {
    'egpd1': ('lkappa',),
    'egpd3': ('ldelta',),
    'egpd4': ('ldelta', 'lkappa'),
    'egpd2direct': ('lkappa1', 'lkappa2', 'logitp'),
}


def multi_indices(n, order):
    out = []
    for combo in itertools.combinations_with_replacement(range(n), order):
        idx = [0] * n
        for c in combo:
            idx[c] += 1
        out.append(tuple(idx))
    return out


def derivative(nll, pars, m):
    e = nll
    for p, k in zip(pars, m):
        for _ in range(k):
            e = sp.diff(e, p)
    return e


def emit_entry(name, entry):
    """C++ helper returning one derivative; entry is 1-based, e.g. (1,2,2)."""
    nll_fn, pars = MODELS[name]
    m = [0] * len(pars)
    for c in entry:
        m[c - 1] += 1
    e = derivative(nll_fn(), pars, m)
    e = e.subs({p: symbols(s, real=True)
                for p, s in zip(pars[2:], PAR_NAMES[name])})
    reps, (red,) = cse([e])
    names = ["y", "lpsi", "xi"] + list(PAR_NAMES[name])
    fn = "%s_d%d_%s" % (name, len(entry), "".join(map(str, entry)))
    print("// GENERATED by derive_egpd_d34.py (entry %s %s) -- do not edit by hand."
          % (name, ",".join(map(str, entry))))
    print("static inline double %s(%s) {"
          % (fn, ", ".join("double %s" % s for s in names)))
    for s, ex in reps:
        print("  const double %s = %s;" % (s, ccode(ex)))
    print("  return %s;" % ccode(red))
    print("}")


def gpd_term():
    """The GPD-density contribution to entry (1,2,2), simplified in z = y/psi."""
    z = symbols('z', positive=True)
    gpd = lpsi + (1 + xi) * log1p(xi * y * exp(-lpsi)) / xi
    d = sp.diff(gpd, lpsi, xi, xi).subs(y, z * exp(lpsi))
    d = sp.factor(sp.simplify(d.rewrite(sp.log)))
    print("d^3 [lpsi + (1 + xi) log1p(xi z)/xi] / (d lpsi d xi^2) =", d)
    print("C++ (inside y * (...)/psi):", ccode(sp.simplify(d / z)).replace('z', '(y/ee1)'))


def cse(exprs):
    syms = sp.numbered_symbols('q')
    return sp.cse(exprs, symbols=syms, optimizations='basic')


CHECK_POINTS = [
    ('egpd1', (3.3, 0.6931471805599453, 0.3, 0.4054651081081644)),
    ('egpd1', (5.1, 0.6931471805599453, -0.2, 0.4054651081081644)),
    ('egpd3', (1.1, 0.6931471805599453, -0.1, -0.6931471805599453)),
    ('egpd3', (3.3, 0.6931471805599453, 0.3, 0.4054651081081644)),
    ('egpd4', (2.5, 0.4, 0.2, -0.3, 0.7)),
    ('egpd2direct', (2.5, 0.4, 0.2, -0.3, 0.7, 0.4)),
]


def check(digits=40):
    for name, vals in CHECK_POINTS:
        nll_fn, pars = MODELS[name]
        nll = nll_fn()
        subs = {y: sp.Rational(repr(vals[0]))}
        for p, v in zip(pars, vals[1:]):
            subs[p] = sp.Rational(repr(v))
        order = multi_indices(len(pars), 3) + multi_indices(len(pars), 4)
        out = []
        for mi in order:
            out.append(sp.N(derivative(nll, pars, mi).subs(subs), digits))
        print("## %s at y=%r eta=%r" % (name, vals[0], vals[1:]))
        print("c(%s)" % ", ".join("%.17g" % float(v) for v in out))


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'check'
    if what == 'check':
        check()
    elif what == 'entry':
        emit_entry(sys.argv[2], tuple(int(s) for s in sys.argv[3].split(',')))
    elif what == 'gpd':
        gpd_term()
    else:
        sys.exit("usage: derive_egpd_d34.py check | gpd | entry <model> i,j,k[,l]")
