from sympy.abc import x
from sympy import ZZ, Poly
from numpy import pad as nppad


class NTRU:
    """ Stores information about the generic NTRU parameters """
    def __init__(self, N, p, q):
        self.N = N
        self.p = p
        self.q = q
        self.R = Poly(x**N - 1, x).set_domain(ZZ)


def asPoly(arr):
    """ Convert a NumPy array to a polynomial """
    return Poly(arr[::-1],x).set_domain(ZZ)


def asArr(poly, N):
    """ Convert polynomial to a NumPy array """
    tmp = poly.all_coeffs()[::-1]
    return nppad(tmp, (0, N - len(tmp)))
